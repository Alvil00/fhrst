#coding:utf-8
import numpy
import h5py
import struct
import logging
import argparse
import enum
import os
import os.path
import collections.abc
import pdb

__author__ = 'Vilenskiy Alexey dep. 50'
log = None

TABLE_OF_ETYPES = {}
SOLID_ELEMENTS = {182, 183, 185, 186, 187}
SOLID_ELEMENTS_2D = {182, 183}

FOUR_TEMPLATE = '{num:6}  "{rname:<11}"{a[0]:>7.1f}{a[1]:>11.2f}{a[2]:>11.2f}{a[3]:>11.2f}{a[4]:>11.2f}       0.00       0.00\n'
SIX_TEMPLATE = '{num:6}  "{rname:<11}"{a[0]:>7.1f}{a[1]:>11.2f}{a[2]:>11.2f}{a[3]:>11.2f}{a[4]:>11.2f}{a[5]:>11.2f}{a[6]:>11.2f}\n'

# Функция считывающая аргументы командной строки
def parse_args():
	main_parser = argparse.ArgumentParser(description='HRST fast analog. Approximately 50 times faster than dep 72 junky version')
	main_parser.add_argument('-l','-log', action='store_true', help='save some log information')
	main_parser.add_argument('-v','-verbose', action='store_true', help='verbose mode')
	
	subparsers = main_parser.add_subparsers(help='sub-commands help', title='subcommands', description='extractors of stress tensor', dest='extractor')
	# Parser from rst file
	fromrst_parser = subparsers.add_parser('from-rst', help='extract stress tensor from rst to hcn binary file')
	
	
	fromrst_parser.add_argument('-f', '-forced', action='store_true', help='forced save each 6 component from stress tensor\
	if that has been disabled like in default case, programm save only 4 component of stress if it has not found any stereo element')
	fromrst_parser.add_argument('-a', '-averaged', action='store_true', help='average stress on material borders')
	
	
	exlcluded_items = fromrst_parser.add_mutually_exclusive_group()
	exlcluded_items.add_argument('--nodefile', type=str, default=None, help='extract only nodes wich included in this file')
	exlcluded_items.add_argument('--elemfile', type=str, default=None, help='extract only nodes wich has been owned by pointed elements wich has been contained in this file. WARNING: BE CAREFUL WITH THIS OPTION!')
	
	fromrst_parser.add_argument('--sm', '--start_moment', default=1, type=int, help='start extraction moment - default set as one')
	fromrst_parser.add_argument('--em', '--end_moment', default=0, type=int, help='end extraction moment - default set as zero so it means that script extract all moments till the last')
	#fromrst_parser.add_argument('--et','--extract_type', choises=('stress', 'elastic_strains', 'plastic_strains', 'thermal_strains', 'creeep_strains'), type=str, default='stress', help='extract neseccey tensor from rst')
	#fromrst_parser.add_argument('--mt','--mat_table', default=None, type=str, help='average results on material borders by group of materials wich contained in pointed file')
	
	fromrst_parser.add_argument('--rstfile', help='name of rst file', metavar='file', required=True)
	fromrst_parser.add_argument('--hcnfile', help='name of hcn file', metavar='file', required=True)
	
	
	# Parser from archive file
	fromhcn_parser = subparsers.add_parser('from-hcn',help='extract stress tensor from hcn binary files to tmp text files')
	
	fromhcn_parser.add_argument('-z', '-zero', action='store_true', help='append zero item to each tmp file head')
	fromhcn_parser.add_argument('-i', '-list', action='store_true', help='create file with list of extracted nodes in each material folder')
	fromhcn_parser.add_argument('-r', '-regtimefile', action='store_true', help='add file wich contain sequence of time moments in each regime')
	
	fromhcn_parser.add_argument('--nodefile', type=str, help = "file wich containing neseccery nodes")
	fromhcn_parser.add_argument('--outdir', type=str, help='dir with final tmp files')
	
	
	fromhcn_parser.add_argument('hcnfiles', type=str, nargs='+', help='binary hcn archive file')
	
	r = main_parser.parse_args()
	return r


def extract_nodes(filename=None, com_str=1, format=(int,)):
	if not filename:
		filename = "NODE.TXT"
	node_set = set()
	# void list to apply format logic to standalone values in line
	c = list()
		# count the number of lines in the file and put it to variable named l
		# it is neseccery only for progress bar
	with open(filename, "r") as node_file:
		for num, i in enumerate(node_file):
			# first two strings in standard NODE.txt file which has created by ansalt program
			# consist from comments
			# then we need to skip them
			if num>=com_str:
				buf_str_list = i.strip().split() # we split the node_str to standalone values 
				if len(buf_str_list)<len(format):
					for q in range(len(format) - len(buf_str_list)):
						buf_str_list.append("0")
				for format_num, format_cast in enumerate(format): # and apply to all of them appropriate cast
					c.append(format_cast(buf_str_list[format_num]))
				node_set.add(c[0]) # append number of node and node group (and maybe some another values included in format argument) in a list:new version
				c = list() # clear list which we appended to returnable list
	return node_set


# Класс описывающий элементы для функции extract_tensor_from_rst
class Element:
	__slots__ = ['num', 'mat', 'etype', 'nlist']
	def __init__(self, num, mat, etype, nlist):
		self.num = num
		self.mat = mat
		self.etype = etype
		self.nlist = list(nlist)
		if TABLE_OF_ETYPES.get(self.etype)[0] in SOLID_ELEMENTS_2D:
			if self.nlist[3]==self.nlist[2]:
				self.nlist.pop()
		elif TABLE_OF_ETYPES.get(self.etype)[0]!=187:
			if self.nlist[2]==self.nlist[3] and self.nlist[4]==self.nlist[5]: #Thetra optinon
				self.nlist = self.nlist[:2]
				self.nlist.append(self.nlist[4])
			elif self.nlist[4]==self.nlist[5]:	#Pyramid option
				self.nlist = self.nlist[:4]
			elif self.nlist[2]==self.nlist[3]:  #Prism option
				self.nlist = self.nlist[:3]
				self.nlist.extend(self.nlist[4:8])
			else:
				pass
		else:
			pass
				
	@classmethod
	def frombytes(cls, num, bseq):
		global TABLE_OF_TYPES
		k = struct.unpack('{}i'.format(30), bseq)
		nnodes = 0
		et = TABLE_OF_ETYPES.get(k[1])
		if et[0] not in SOLID_ELEMENTS:
			return None
		return cls(num, k[0], k[1], k[10:(10  + et[1])])


# Класс определяющий структуру узла для метода extract_tensor_from_rst
class Node:
	__slots__ = ['num', 'tmtable', '_num_of_parents']
	def __init__(self, num, nsets, base):
		self.num=num
		self.tmtable = numpy.zeros((nsets, base), dtype=numpy.single)
		self._num_of_parents = 0
	
	def append_time_moment(self, timemoment, stress_tensor):
		self.tmtable[timemoment] += stress_tensor
		
	def norm(self):
		self.tmtable = self.tmtable.dot(1 / self._num_of_parents)
	
	def num_of_parents_addone(self):
		self._num_of_parents += 1

	
class MatZone:
	def __init__(self, matnum):
		self.nlist = {}
		self.mat = matnum

		
class ExcludedType(enum.Enum):
	NON = 0
	NODE = 1
	ELEMENT = 2


# Возвращает список позиций элементов sub_list в листе base_list и сам элемент base_list
def find_index_list(sub_list, base_list):
	result = []
	for i in sub_list:
		result.append((base_list.index(i), i))
	result.sort(key=lambda a:a[0])
	return result


# Экспорт итоговой стрктуры работы функции extract_tensor_from_rst в файл hdf		
def export_to_mat_zone_bin_file(filename, mat_struct):
	f = h5py.File(filename, mode='w')
	for zone in mat_struct.keys():
		nrf = f.create_group(str(zone))
		for node in mat_struct[zone].nlist.keys():
			nrf.create_dataset(name=str(node), data=mat_struct[zone].nlist[node].tmtable)
	f.flush()
	f.close()


# Функция которая возвращает из опредленного набора step-substep-n_last_iter только псоледние
# сошедшиеся substeps c номером set'a
def last_substep_generator(ls, p_steps):
	k = 0
	for n, i in enumerate(ls):
		if i[1] <= k:
			yield n - 1, p_steps[n-1]
		k = i[1]
	yield n, p_steps[n]


# Функция извлечения тензора из фалов типа rst и опционально запись ее в файл типа hdf
# rstname - имя входного файла разрешения *.rst
# hсnname - имя выходного файла (если None - то возвращает структуру)
# base - количество независмых ненулевых компонент в извлекаемом тензоре (4 - при духмерном случае, 6 - при произвольном случае)
# excluded_set - исключенное множество элементов или узлов для которых будет извлечен тензор
# excluded_type - тип исключенного множества (элементы или узлы)
# average - флаг отвечающий за усреднение на границах раздела материалов
# verbose - выводит больше информации
# start_moment - начальный loadstep извлечения
# end_moment - конечный loadstep извлечения (если 0 - извлекает до посленего loadstep).
def extract_tensor_from_rst(rstname='file.rst', hcnname=None, base=None, excluded_set=None, excluded_type:ExcludedType=ExcludedType.NON, average=False, verbose=False, start_moment=1, end_moment=0):
	global log
	with open(rstname, 'rb', buffering=0) as f:
		# Этот участок для чтения частей заголовка
		f.seek(12)
		d = f.read(4)	
		file_format = struct.unpack('i', d)[0]
		
		f.seek(28, 1)
		d = f.read(4)
		ansys_version = float(struct.unpack('4s', d)[0][::-1])
		
		f.seek(8 + 26 * 4, 0)
		d = f.read(4)
		max_file_length = struct.unpack('i', d)[0]
		
		f.seek(424)
		d = f.read(36)
		maxn, nnod, resmax, numdof, maxe, nelm, kan, nsets, ptrend = struct.unpack('9i', d)
		
		d = f.read(24)
		ptr_data_steps, ptr_time_val, ptr_ls, ptr_elems, ptr_nodes, ptr_geom = struct.unpack('6i', d)
		
		# Список узлов
		f.seek(ptr_nodes * 4 + 8)
		d = f.read(4 * nnod)
		nodes = list(struct.unpack('{}i'.format(nnod), d))
		nodes.sort()

		# Список элементов
		f.seek(ptr_elems * 4 + 8)
		d = f.read(4 * nelm)
		elements = list(struct.unpack('{}i'.format(nelm), d))
		
		# Список шагов
		f.seek(ptr_data_steps * 4 + 8)
		d = f.read(4 * nsets)
		ptr_steps = list(struct.unpack('{}I'.format(nsets), d))
		max_int_size_step = ptr_steps[0]
		
		# Поиск шага при котором указатель на него становится длиным (длинный указатель вычиляется как lptr = 2**34 + ptr)
		for n, i in enumerate(ptr_steps):
			if max_int_size_step>i:
				max_int_size_step = n + 1
				break
			max_int_size_step = i
		
		# Таблица шагов, подшагов и итераций - нужна для определения последних подшагов
		f.seek(ptr_ls * 4 + 8)
		d = f.read(12 * nsets)
		ls_array = numpy.frombuffer(d, dtype=numpy.int32)
		ls_array = numpy.reshape(ls_array, (nsets, 3))
		ls_final_list = list(last_substep_generator(ls_array, ptr_steps))
		nfsets = len(ls_final_list)
	
		# Доступ к геометрии
		f.seek((ptr_geom) * 4 + 8)
		d = f.read(180)
		k = struct.unpack('45i', d)
		n_etypes = k[1]
		ptr_etypes = k[20] + k[21]
		size = k[18]
		
		# Если работает усреднение на границах материалов (аргумент функции average) просто представляем все элементы как один материал
		if not average:
			n_mat = k[13]
		else:
			n_mat = 1
		
		# Находим длинный указатель на таблицу элементов
		ptr_eid = k[28] + k[29]
		
		# Логирование информации
		if log:
			log.info('ff\t {}'.format(file_format))
			log.info('matnum\t {}'.format(n_mat))
			log.info('aver\t {}'.format(ansys_version))
			log.info('flength\t {}'.format(max_file_length))
			log.info('nsets\t {}'.format(nsets))
			log.info('nnodes\t {}'.format(nnod))
			log.info('nelms\t {}'.format(nelm))
			log.info('table of sets:')
			previous_set = (0, 0, 0)
			# Алгоритм вывода таблици set'ов
			if len(ls_array)>1:
				for current_set in ls_array:
					if current_set[1]<=previous_set[1]:
						pulse = '├─>'
					else:
						pulse = '│  '
					if previous_set[0]:
						log.info("{pulse}{o[0]:5}{o[1]:5}{o[2]:5}".format(o=previous_set, pulse=pulse))
					previous_set = current_set
			else:
				current_set = ls_array[0]
			log.info((u"└─>{o[0]:5}{o[1]:5}{o[2]:5}".format(o=current_set)))
		
		# Доступ к типам элементов
		f.seek(ptr_etypes * 4 + 8)
		d = f.read(4 * n_etypes)
		etypes_ptr_list = list(struct.unpack('{}i'.format(n_etypes), d))
		for i in etypes_ptr_list:
			f.seek((i + ptr_etypes) * 4 + 8)
			d = f.read(8)
			etn = struct.unpack('ii', d)
			f.seek(364, 1)
			d = f.read(4)
			nnodes_per_element_ws = struct.unpack('i', d)
			TABLE_OF_ETYPES[etn[0]] = (etn[1], nnodes_per_element_ws[0])
		
		# Назначаем количество ненулевых элементов в извлекаемом тензоре при отсутствии назначенного
		# Если значение base при передачи в аргументы функции находится не в допускаемом множестве значений (4,6) то пытаемся сами назначить значение base
		if base is None or base not in (4, 6):
			# В случае если все типы элементов плоские(2 степени свободы на узел) выбираем 4 иначе 6
			if all(map(lambda a: a[0] in SOLID_ELEMENTS_2D and a[0] in SOLID_ELEMENTS or a[0] not in SOLID_ELEMENTS, TABLE_OF_ETYPES.values())):
				base = 4
			else:
				base = 6
		
		base_plus_temp = base + 1
		
		if excluded_type == ExcludedType.ELEMENT:
			elements_walk = find_index_list(excluded_set, elements)
	
			
		# Строим таблицу элементов
		f.seek(ptr_eid * 4 + 8)
		d = f.read(8 * nelm)
		ptr_elem_table = struct.unpack('{}q'.format(nelm), d)
		enod_table = {}
		
		# Определяем element walk - способ прохождения по всем исключенным элементам для разных типов excluded set
		if excluded_type == ExcludedType.NON or excluded_type == ExcludedType.NODE:
			for enum, etable_ptr in zip(elements, ptr_elem_table):
				f.seek((ptr_eid + 2 + etable_ptr) * 4)
				d = f.read(120)
				_e = Element.frombytes(enum, d)
				if _e is not None:
					if average:
						_e.mat = 1
					enod_table[enum] = _e
		elif excluded_type==ExcludedType.ELEMENT:
			for enumintable, enum in elements_walk:
				f.seek((ptr_eid + 2 + ptr_elem_table[enumintable]) * 4)
				d = f.read(120)
				_e = Element.frombytes(enum, d)
				if _e is not None:
					enod_table[enum] = _e
		
		if excluded_type==ExcludedType.NODE:
			elements_walk = []
			for enumintable, enum in enumerate(elements):
				cur_elem = enod_table.get(enum, None)
				if cur_elem is not None and any(map(lambda x:x in excluded_set, cur_elem.nlist)):
					elements_walk.append((enumintable, enum))		

				
		# Назначаем начальный извлекаемый момент времени и конечный извлекаемый момент времени
		if start_moment==1 and end_moment==0:
			time_moment_check_ranage = None
		else:
			if end_moment==0:
				end_moment = nfsets
			if end_moment<=start_moment:
				raise IndexError("Wrong end moment: it is lesser then start moment")
			if end_moment>nfsets:
				raise IndexError("Wrong end moment: it is bigger than file last moment")
			if start_moment<=0:
				raise IndexError("Wrong start moment: it is lesser than zero or equal to then")
			time_moment_check_ranage = range(start_moment, end_moment + 1)
			nfsets = len(time_moment_check_ranage)
		
		# Итоговая результирующая структура которая будет хранить финальные результаты
		result_structure = {}
		
		# Заполняем ее материальными зонами
		for i in range(1, n_mat + 1):
			result_structure[i] = MatZone(i)
			
		# При привышении указателя размером 2**34 к нему надо добавлять это число
		add_size = 0
		_add_size = 17179869184 #2**34
		
		# Переменная curnum используется чтобы в случае несоответствии количества режимов при назначенных startmoment и endmoment 
		# режимы извекались под правильным номером
		curnum = 0
		
		# Основной цикл который извлекает тензор и температуры на каждом шаге
		for step_num, (set_num, step_ptr) in enumerate(ls_final_list, 1):
			if set_num>=max_int_size_step - 1:
				add_size = _add_size
			f.seek(step_ptr * 4 + 8 + add_size)  # Переходим к заголовку режима
			d = f.read(48)			# Считываем заголовок	
			s_header = struct.unpack('12i', d)
			extracted_step_num = s_header[4]		# Текущий номер режима извлеченный из таблицы
			esol_ptr = s_header[11]	# Указатель на element solution items (указатели измеряются от начала текущего режима)
			if extracted_step_num==step_num:
				if time_moment_check_ranage is not None and extracted_step_num not in time_moment_check_ranage:
					continue
				if verbose:
					print("timemoment %i" % extracted_step_num)
				esol_ptr_int = (esol_ptr - 12) * 4
				f.seek(esol_ptr_int, 1) # Перемещаемся к таблицы в которой содержаться указатели на решение в каждом из элеметов
				# По всем элементам
				d = f.read(8 * nelm) # Считываем ленту элементов (ряд указателей на реузльтаты)
				f.seek(-8 * nelm, 1)
				h = struct.unpack('{}q'.format(nelm), d)
				if excluded_type == ExcludedType.NON:
					elements_walk = enumerate(elements)
				for progress, elementnum in elements_walk:
					element_instance = enod_table.get(elementnum, None)
					if element_instance is None:
						continue
					cur_nlist = result_structure[element_instance.mat].nlist
					ptr_to_solution_info = h[progress]
					f.seek(ptr_to_solution_info * 4, 1) # Переходим по указателю к информации о элементе
					d = f.read(12) # Считываем заголовки информации о элементе
					ptr_to_nodal_stresses = struct.unpack('3i', d)
					ptr_to_nodal_stresses = ptr_to_nodal_stresses[2]
					f.seek(52, 1) #16*4-3*4 где 3 - количество int в заголовке который мы считали строками выше
					d = f.read(4)
					ptr_to_temp = struct.unpack('i', d)[0]
					f.seek(ptr_to_nodal_stresses * 4 - 68, 1) #Переход к полям напряжений для узла
					stresses_ = {}
					# Вынимаем необходимый тензор для всех узлов элемента
					for corner_node_num in element_instance.nlist:
						corner_node = cur_nlist.get(corner_node_num, None)
						if corner_node is None:	
							cur_nlist[corner_node_num] = Node(corner_node_num, nfsets, base_plus_temp)
							corner_node = cur_nlist.get(corner_node_num)
						if step_num==1:
							corner_node.num_of_parents_addone()
						d = f.read(24)
						stresses_[corner_node_num] = d
					f.seek((ptr_to_temp - ptr_to_nodal_stresses - 6 * len(element_instance.nlist)) * 4, 1) #Переход к полю температуры для узла
					# Вынимаем температуру и записываем значение для каждого узла
					for corner_node_num in element_instance.nlist:
						corner_node = cur_nlist.get(corner_node_num)
						d = f.read(4)
						# Получается здесь мы складываем(concatenate) две байтовые строки а затем восстанавливаем их прямо в numpy archive с типом numpy.single (32 bit float)
						stresses_[corner_node_num] = numpy.frombuffer(d + stresses_[corner_node_num], dtype=numpy.single, count=base_plus_temp)
						corner_node.append_time_moment(curnum, stresses_[corner_node_num])
					f.seek((-ptr_to_temp - len(element_instance.nlist) - ptr_to_solution_info) * 4, 1) # Возврат к системе координат
				curnum += 1
			else:
				raise IndexError("Wrong time moment {} at number {}.Programm will be terminated".format(extracted_step_num, step_num))
				
		# Усредняем знчаения напряжений для элементов
		if excluded_type != ExcludedType.NODE:
			for zone in result_structure.keys():
				for node in result_structure[zone].nlist.values():
					node.norm()
		else:
			# В случае если тип услючения у нас стоит по узлам, мы все равно пробегали по всем узлам что принадлежат их родительским элементам.
			# Таким образом необходимо отсеять узлы которые не входили в excluded set
			for zone in result_structure.keys():
				lost_nodes = set(result_structure[zone].nlist.keys()).difference(excluded_set)
				for node in lost_nodes:
					result_structure[zone].nlist.pop(node)
				for node in result_structure[zone].nlist.values():
					if node.num in excluded_set:
						node.norm()
		
		
		# Экспорт в файл или вернуть как структуру в случае если outname is None
		if hcnname is not None:
			export_to_mat_zone_bin_file(hcnname, result_structure)
			return None
		else:
			return result_structure


# Функция записывает файл в котором есть режимы с последовательностью моментов вермени
def create_regtime_file(file, sequences):
	if all(map(lambda a:isinstance(a, collections.abc.Iterable), sequences)):
		with open(file, 'w') as f:
			for n, i in enumerate(sequences, 1):
				f.write('r{n}_1 "r\t\t" {s}\n'.format(n=n, s=str(list(i))[1:-1].replace(',', '')))
	

# Функция извлекает тензор напряжения из фала типа hcn(hdf) в файлы типа tmp
# filelist - список имен(или путей) hcn(hdf) файлов
# nodefile - файл с номерами узлов которые необходимо извлеч(первая строка в этом файле комментарий). Елси значение None - извлекает все узлы.
# outdir - папка в которую извлекаются узлы
# verbose - выводить в консоль больше информации
# zeros - флаг который регулирует добавление нулевого момента времени к файлу tmp
# listfile - флаг который создает в каждой папке материалов файл NODE.TXT со списком узлов в данной папке
# regtime - флаг регулирующий содание файла реимов (файл n.his в папке запуска скрипта)
def extract_tensor_from_hcn(filelist, nodefile=None, outdir=None, verbose=False, zeros=False, listfile=False, regtime=False):
	global log
	# Словарь из элемнтов:ключ - номер материала, значение - множество узлов в этом материале
	mat_node_set = {}

	if nodefile is None:
		nodeset = None
	else:
		nodeset = extract_nodes(nodefile)
		if verbose:
			print("Length of nodefile:{}".format(len(nodeset)))
	if outdir is None:
		outdir = os.curdir
	curnum = 1
	
	temp = None
	
	if regtime:
		regtimelist = []
	
	# Для всех файлов в списке файлов hcn
	for file in filelist:
		wmode = 'a' if curnum > 1 else 'w'
		regime_name = os.path.basename(file).split('.')[0][:11]
		if verbose:
			print("Try to unpack {}".format(os.path.basename(file)))
		# Пытаемся открыть файл на чтение
		with h5py.File(file, mode='r') as f:
			# Для каждого материала в файле
			for mat in f.keys():
				# На самом первом шаге для всех файлов создаем необходимые папки для извлекаемой информации
				if curnum==1:
					# Определяем множество mat_node_set как ПЕРЕСЕЧЕНИЕ множества узлов в множестве nodeset и множества узлов в файле
					if nodeset is not None:
						mat_node_set[mat] = nodeset.intersection(set(map(int, f[mat].keys())))
					# Если nodeset отсутвсвует просто берем набор узлов в файле как множество mat_node_set
					else:
						mat_node_set[mat] = set(map(int, f[mat].keys()))
					if len(mat_node_set[mat]):
						try:
							os.mkdir('{o}{sep}MATERIAL_{m}'.format(m=mat, sep=os.path.sep, o=outdir))
						except FileExistsError:
							pass
				# Для каждого узла который входит в множество mat_node_set
				for node in mat_node_set[mat]:
					# Записываем в файл (файл открывается в моде 'w' при curnum==1 и в моде 'a' во всех остальных случаях)
					of = open('{o}{sep}MATERIAL_{m}{sep}{n}n.tmp'.format(n=node, m=mat, sep=os.path.sep, o=outdir), mode=wmode)
					k = numpy.array(f[mat][str(node)])
					# При флаге zeros добавляем в файл строку с нулевым моментом времени
					if curnum == 1 and zeros:
						of.write(FOUR_TEMPLATE.format(a = (20.0, 0.0, 0.0, 0.0, 0.0), num=0, rname="NULL"))
					# Для каждого момента времени в узле в файле
					for n,i in enumerate(k, curnum):
						if temp is None and n==1:
							if k.shape[1]==5:
								temp = FOUR_TEMPLATE
							else:
								temp = SIX_TEMPLATE
						# Записываем строку с тензором
						of.write(temp.format(a=i, num=n, rname=regime_name))
					of.close()
		if regtime:
			regtimelist.append(range(curnum, n + 1))
		# Принимаем curnum n+1 для того чтобы момент времени из следующего файла пошел под сквозной нумерацией 
		curnum = n + 1
	if listfile:
		for mat in mat_node_set.keys():
			if len(mat_node_set[mat]):
				with open("{o}{sep}MATERIAL_{m}{sep}ZONA1.TXT".format(m=mat, sep=os.path.sep, o=outdir), mode='w') as zf:
					zf.write("%i\n" % len(mat_node_set[mat]))
					for i in sorted(mat_node_set[mat]):
						zf.write("%i\n" % i)
	if regtime:
		create_regtime_file('n.his', regtimelist)


def main():
	global log
	args = parse_args()
	if args.l and args.extractor is not None:# Включено логирование
		fhand = logging.FileHandler(filename='fhrstlog.out', mode='w', encoding='utf-8')
		fhand.setLevel(level=0)
		formstr = "%(message)s"
		log = logging.Logger("mainlog", level=0)
		log.addHandler(fhand)
		log.info("=-----LOG START-----=")
		if args.v and log:
			log.info('additional information will be output in console.')
	if args.extractor=='from-rst':
		base = None
		excluded_type = ExcludedType.NON
		if args.elemfile is not None:
			excluded = extract_nodes(args.elemfile)
			excluded_type = ExcludedType.ELEMENT
			if log:
				log.info('extract nodes wich owned by elements wich lie in {} file'.foramt(args.elemfile))
		elif args.nodefile is not None:
			excluded = extract_nodes(args.nodefile)
			excluded_type = ExcludedType.NODE
			if log:
				log.info('extract nodes wich contain in {} file'.foramt(args.nodefile))
		else:
			excluded = None
			if log:
				log.info('no exclusion for any item')
		if args.f:
			base = 6
		if args.a and log:
			log.info('nodes on border between materials will be averaged by value')
		extract_tensor_from_rst(args.rstfile, args.hcnfile, base, excluded, excluded_type, args.a, args.v, args.sm, args.em)
	elif args.extractor=='from-hcn':
		extract_tensor_from_hcn(args.hcnfiles, args.nodefile, args.outdir, args.v, args.z, args.i, args.r)
	else:
		print('print -h or --help key to show help message.')
	
if __name__=="__main__":
	main()
