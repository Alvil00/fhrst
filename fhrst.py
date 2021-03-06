#coding:utf-8
import numpy
import h5py
import struct
import logging
import argparse
import enum
import os
import os.path
import sys
import collections
import collections.abc


try:
	import unite
	GLOBAL_UNITE_LACK_MESSAGE = ''
except ImportError:
	unite = None
	GLOBAL_UNITE_LACK_MESSAGE = '!WARNING! This routine is not available because python cant find out unite module'
	


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
	
	subparsers = main_parser.add_subparsers(help='sub-commands help', title='subcommands', description='extractors of stress tensor and another useful features', dest='extractor')
	
	# Parser from rst file
	fromrst_parser = subparsers.add_parser('from-rst', help='extract stress tensor from rst to hcn binary file')
	
	fromrst_parser.add_argument('-f', '-forced', action='store_true', help='forced save each 6 component from stress tensor\
	if that has been disabled like in default case, programm save only 4 component of stress if it has not found any stereo element')
	fromrst_parser.add_argument('-a', '-averaged', action='store_true', help='average stress on material borders')
	
	
	exlcluded_items = fromrst_parser.add_mutually_exclusive_group()
	exlcluded_items.add_argument('--nodefile', type=str, default=None, help='extract only nodes included in this file')
	exlcluded_items.add_argument('--elemfile', type=str, default=None, help='extract only nodes have been belonged pointed elements have been contained in this file. WARNING: BE CAREFUL WITH THIS OPTION!')
	
	fromrst_parser.add_argument('--sm', '--start-moment', default=1, type=int, help='start extraction moment - default set as one')
	fromrst_parser.add_argument('--em', '--end-moment', default=0, type=int, help='end extraction moment - default set as zero so it means that script extract all moments till the last')
	#fromrst_parser.add_argument('--et','--extract-type', choices=('stress', 'elastic_strains', 'plastic_strains', 'thermal_strains', 'creeep_strains'), type=str, default='stress', help='extract neseccey tensor from rst')
	#fromrst_parser.add_argument('--mt','--mat-table', default=None, type=str, help='average results on material borders by group of materials which contained in pointed file')
	
	fromrst_parser.add_argument('--rstfile', help='name of rst file', metavar='file', required=True)
	fromrst_parser.add_argument('--hcnfile', help='name of hcn file', metavar='file', required=True)
	
	
	# Parser from archive file
	fromhcn_parser = subparsers.add_parser('from-hcn',help='extract stress tensor from hcn binary files to tmp text files')
	
	zero = fromhcn_parser.add_mutually_exclusive_group()
	zero.add_argument('-z', '-zero', action='store_true', help='append zero item to each tmp file head')
	zero.add_argument('-t', '-tail', action='store_true', help='append zero item to each tmp file tail')
	
	fromhcn_parser.add_argument('-i', '-list', action='store_true', help='create file with list of extracted nodes in each material folder')
	fromhcn_parser.add_argument('-r', '-regtimefile', action='store_true', help='add file which contain sequence of time moments in each regime')
	
	compress_group = fromhcn_parser.add_mutually_exclusive_group()
	compress_group.add_argument('-c', '-check-acompress', action='store_true', help='check for every node - if all calculation statement stresses lesser or equal (node always in absolute compression) then zero then script hasnt create tmp file')
	compress_group.add_argument('-C', '-check-acompress-macro', action='store_true', help="write only the macro which contain zero values of nodes which stay in absolute compression, don't create *.tmp files, you can't use this option with -z, -t, -i, --outdir, -r option", dest='macro_compress')
	fromhcn_parser.add_argument('--set-check-acompress-method', type=str, default='all', choices=('normal', 'all'), help='set which determinate what components we must evaluate for checking, all or only normal', dest='ca_method')
	fromhcn_parser.add_argument('--set-check-acompress-tolerance', type=float, default=0.0 , help='set value of overtop for stress in absolute compress checking', dest='ca_tol', metavar='TOLERANCE')
	check_acompress_method_dict = {'all': CheckAbsoluteCompressMethod.ALL, 'normal':CheckAbsoluteCompressMethod.NORMAL}
	
	fromhcn_parser.add_argument('--nodefile', type=str, help = "file which containing neseccery nodes")
	fromhcn_parser.add_argument('--outdir', type=str, help='dir with final tmp files')
	
	
	fromhcn_parser.add_argument('hcnfiles', type=str, nargs='+', help='binary hcn archive file')
	
	# Loader to rst file
	if unite:
		torst_help_message = 'load accumulated damage into rst file into sx item of nodes results'
	else:
		torst_help_message = GLOBAL_UNITE_LACK_MESSAGE
		
	torst_loader = subparsers.add_parser('to-rst', help=torst_help_message)
	
	torst_loader.add_argument('--rstfile', type=str, help='name of rst file with neseccery geometry', metavar='file', required=True)
	torst_loader.add_argument('--damfile', type=str, help='name of new accumulated damage file', metavar='file', required=True)
	
	torst_loader.add_argument('--mdf', '--mat-damage-file', action='append', dest='mdf', nargs=2,
								help='append material damage macro file. ness 2 args - {num of material} {name of file}',
								metavar=("material_num", "damage_file"), required=True)
	
	#torst_loader.add_argument('-a', '-averaged', action='store_true', help='if this flag is active chose material with try to apply damage from material number 0')
	r = main_parser.parse_args()
	if r.extractor == 'from-hcn':
		r.ca_method = check_acompress_method_dict[r.ca_method]
	ATTRIBUTE_ERROR_TEMPLATE = sys.argv[0] + ": error: you can't use {forrbiden_options} option with option {usable_option}"
	try:
		if r.extractor == 'from-hcn' and r.macro_compress:
			if unite is None:
				raise ImportError
			if r.z:
				raise AttributeError(ATTRIBUTE_ERROR_TEMPLATE.format(forbidden_option='-z or -zero', usable_option='-C'))
			if r.t:
				raise AttributeError(ATTRIBUTE_ERROR_TEMPLATE.format(forbidden_option='-t or -tail', usable_option='-C'))
			if r.i:
				raise AttributeError(ATTRIBUTE_ERROR_TEMPLATE.format(forbidden_option='-i or -list', usable_option='-C'))
			if r.r:
				raise AttributeError(ATTRIBUTE_ERROR_TEMPLATE.format(forbidden_option='-r or -regtimefile', usable_option='-C'))
			if r.outdir:
				raise AttributeError(ATTRIBUTE_ERROR_TEMPLATE.format(forbidden_option='--outdir', usable_option='-C'))
	except AttributeError as ae:
		print(ae.args[0])
		assert True
	except ImportError:
		print(GLOBAL_UNITE_LACK_MESSAGE)
		print('To use -C option put unite module in neseccery directory')
	
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
		elif TABLE_OF_ETYPES.get(self.etype)[0]!=187 and TABLE_OF_ETYPES.get(self.etype)[0] in SOLID_ELEMENTS:
			if self.nlist[2]==self.nlist[3] and self.nlist[4]==self.nlist[5]: #Thetra optinon
				self.nlist = self.nlist[:2] + [self.nlist[4],]
			elif self.nlist[4]==self.nlist[5]:	#Pyramid option
				self.nlist = self.nlist[:4]
			elif self.nlist[2]==self.nlist[3]:  #Prism option
				self.nlist = self.nlist[:3] + self.nlist[4:7]
			else:
				pass
		else:
			pass
	
	# Метод для инициализации элементов из последовательности бинарных данных по 30 int
	@classmethod
	def frombytes(cls, num, bseq):
		global TABLE_OF_ETYPES
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

# Класс в котором хранится номер материала и все принадлежащие к нему узлы
class MatZone:
	def __init__(self, matnum):
		self.nlist = {}
		self.mat = matnum

		
class IncludedType(enum.Enum):
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
	k = (0, 0, 0)
	for n, i in enumerate(ls):
		if i[1] <= k[1] or (i[0] > k[0] and n != 0):
			yield n - 1, p_steps[n-1]
		k = i
	yield n, p_steps[n]


# Функция извлечения тензора из фалов типа rst и опционально запись ее в файл типа hdf
# rstname - имя входного файла разрешения *.rst
# hсnname - имя выходного файла (если None - то возвращает структуру)
# base - количество независмых ненулевых компонент в извлекаемом тензоре (4 - при духмерном случае, 6 - при произвольном случае)
# included_set - включенное множество элементов или узлов для которых будет извлечен тензор
# included_type - тип включенного множества (элементы или узлы)
# average - флаг отвечающий за усреднение на границах раздела материалов
# verbose - выводит больше информации
# start_moment - начальный loadstep извлечения
# end_moment - конечный loadstep извлечения (если 0 - извлекает до посленего loadstep).
def extract_tensor_from_rst(rstname='file.rst', hcnname=None, base=None, included_set=None, included_type:IncludedType=IncludedType.NON, average=False, verbose=False, start_moment=1, end_moment=0):
	global log
	global TABLE_OF_ETYPES
	if verbose:
		print("From file {}".format(rstname))
	with open(rstname, 'rb', buffering=0) as f:
		# Чтение основного бинарного заголовка ANSYS
		if verbose:
			print('Start to read header data...')
		info_header = RTable(f)
		file_format = info_header[1]
		ansys_version = float(struct.unpack('4s', struct.pack('i', info_header[9]))[0][::-1])
		max_file_length = info_header[26]
		
		
		# Чтение заголовка rst
		rst_header = RTable(f)
		maxn, nnod, resmax, numdof, maxe, nelm, kan, nsets, ptrend = rst_header[1:10]
		ptr_data_steps, ptr_time_val, ptr_ls, ptr_elems, ptr_nodes, ptr_geom = rst_header[10:16]

		# Список элементов
		f.seek(ptr_elems * 4)
		elements = RTable(f)
		
		# Список шагов
		f.seek(ptr_data_steps * 4)
		ptr_steps = RTable(f)
		max_int_size_step = ptr_steps[0]
		
		# Поиск шага при котором указатель на него становится длиным (длинный указатель вычиляется как lptr = 2**34 + ptr)
		pstep = ptr_steps[0]
		miszs = []
		for n, i in enumerate(ptr_steps[1:],1):
			if pstep > i:
				miszs.append(n + 1)
			pstep = i
		max_int_size_step = iter(miszs)
		
		# Таблица шагов, подшагов и итераций - нужна для определения последних подшагов
		f.seek(ptr_ls * 4)
		ls_array = numpy.array(RTable(f)[0:nsets * 3])
		ls_array = numpy.reshape(ls_array, (nsets, 3))
		ls_final_list = list(last_substep_generator(ls_array, ptr_steps))
		nfsets = len(ls_final_list)
	
		# Доступ к геометрии
		f.seek(ptr_geom * 4)
		geom_header = RTable(f)
		n_etypes = geom_header[1]
		size = geom_header[18]
		ptr_etypes = geom_header[20] + geom_header[21]
		
		# Если работает усреднение на границах материалов (аргумент функции average) просто представляем все элементы как один материал
		if not average:
			n_mat = geom_header[13]
		else:
			n_mat = 1
		
		# Находим длинный указатель на таблицу элементов
		ptr_eid = geom_header[28] + geom_header[29]
		
		# Логгирование информации
		if log:
			if verbose:
				print('Write some log information...')
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
				for n_cs, current_set in enumerate(ls_array):
					if current_set[1] <= previous_set[1] or (previous_set[0] < current_set[0] and n_cs != 0):
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
		f.seek(ptr_etypes * 4)
		etypes_ptr_list = RTable(f)
		for ptr_local in etypes_ptr_list:
			f.seek((ptr_local + ptr_etypes) * 4)
			current_etype = RTable(f)
			etn = current_etype[0:2]
			nnodes_per_element_ws = current_etype[93]
			TABLE_OF_ETYPES[etn[0]] = (etn[1], nnodes_per_element_ws)
		
		# Назначаем количество ненулевых элементов в извлекаемом тензоре при отсутствии назначенного
		# Если значение base при передачи в аргументы функции находится не в допускаемом множестве значений (4,6) то пытаемся сами назначить значение base
		if base is None or base not in (4, 6):
			# В случае если все типы элементов плоские(2 степени свободы на узел) выбираем 4 иначе 6
			if all(map(lambda a: a[0] in SOLID_ELEMENTS_2D and a[0] in SOLID_ELEMENTS or a[0] not in SOLID_ELEMENTS, TABLE_OF_ETYPES.values())):
				base = 4
			else:
				base = 6
		
		# Необходимо чтобы не вычислять это значение каждый раз в цикле
		base_plus_temp = base + 1
		
		if included_type == IncludedType.ELEMENT:
			elements_walk = find_index_list(included_set, elements)
	
		if verbose:
			print('Build element table...')
		# Строим таблицу элементов
		f.seek(ptr_eid * 4)
		ptr_elem_table = RTable(f, typ='q')
		enod_table = {}
		
		# Определяем element walk - способ прохождения по всем включенным элементам для разных типов included set
		if included_type == IncludedType.NON or included_type == IncludedType.NODE:
			for enum, etable_ptr in zip(elements, ptr_elem_table):
				f.seek((ptr_eid + 2 + etable_ptr) * 4)
				d = f.read(120)
				_e = Element.frombytes(enum, d)
				if _e is not None:
					if average:
						_e.mat = 1
					enod_table[enum] = _e
		elif included_type==IncludedType.ELEMENT:
			for enumintable, enum in elements_walk:
				f.seek((ptr_eid + 2 + ptr_elem_table[enumintable]) * 4)
				d = f.read(120)
				_e = Element.frombytes(enum, d)
				if _e is not None:
					enod_table[enum] = _e
		
		if included_type==IncludedType.NODE:
			elements_walk = []
			for enumintable, enum in enumerate(elements):
				cur_elem = enod_table.get(enum, None)
				if cur_elem is not None and any(map(lambda x:x in included_set, cur_elem.nlist)):
					elements_walk.append((enumintable, enum))		

				
		# Назначаем начальный извлекаемый момент времени и конечный извлекаемый момент времени
		if start_moment==1 and end_moment==0:
			time_moment_check_ranage = None
			end_moment = nfsets
		else:
			if end_moment==0:
				end_moment = nfsets
			if end_moment < start_moment:
				raise IndexError("Wrong end moment: it is lesser then start moment")
			if end_moment > nfsets:
				raise IndexError("Wrong end moment: it is bigger than file last moment")
			if start_moment <= 0:
				raise IndexError("Wrong start moment: it is lesser than zero or equal to then")
			time_moment_check_ranage = range(start_moment, end_moment + 1)
			nfsets = len(time_moment_check_ranage)
		if verbose:
			print("Start moment:\t{:>7}".format(start_moment))
			print("End moment:\t{:>7}".format(end_moment))
		# Итоговая результирующая структура которая будет хранить финальные результаты
		result_structure = {}
		
		# Заполняем ее материальными зонами
		for i in range(1, n_mat + 1):
			result_structure[i] = MatZone(i)
			
		# При привышении указателя размером 2**34 к нему надо добавлять это число
		add_size = 0
		_add_size = 17179869184 #2**34
		if verbose:
			print("Start the main read cycle...")
		# Переменная curnum используется чтобы в случае несоответствии количества режимов при назначенных startmoment и endmoment 
		# режимы извекались под правильным номером
		curnum = 0
		current_max_int_size_step = next(max_int_size_step)
		# Основной цикл который извлекает тензор и температуры на каждом шаге
		for step_num, (set_num, step_ptr) in enumerate(ls_final_list, 1):
			if set_num >= current_max_int_size_step - 1:
				add_size += _add_size
				current_max_int_size_step = next(max_int_size_step)
			f.seek(step_ptr * 4 + 8 + add_size)  # Переходим к заголовку режима
			d = f.read(48)			# Считываем заголовок	
			s_header = struct.unpack('12i', d)
			extracted_step_num = s_header[4]		# Текущий номер режима извлеченный из таблицы
			esol_ptr = s_header[11]	# Указатель на element solution items (указатели измеряются от начала текущего режима)
			if extracted_step_num == step_num:
				if time_moment_check_ranage is not None and extracted_step_num not in time_moment_check_ranage:
					continue
				if verbose:
					print("Current time moment read {}/{}\t".format(extracted_step_num, end_moment),end='\r')
				esol_ptr_int = (esol_ptr - 12) * 4
				f.seek(esol_ptr_int, 1) # Перемещаемся к таблицы в которой содержаться указатели на решение в каждом из элеметов
				# По всем элементам
				d = f.read(8 * nelm) # Считываем ленту элементов (ряд указателей на реузльтаты)
				f.seek(-8 * nelm, 1)
				h = struct.unpack('{}q'.format(nelm), d)
				if included_type == IncludedType.NON:
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
						if step_num == start_moment:
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
		if included_type != IncludedType.NODE:
			for zone in result_structure.keys():
				for node in result_structure[zone].nlist.values():
					node.norm()
		else:
			# В случае если тип услючения у нас стоит по узлам, мы все равно пробегали по всем узлам что принадлежат их родительским элементам.
			# Таким образом необходимо отсеять узлы которые не входили в included set
			for zone in result_structure.keys():
				lost_nodes = set(result_structure[zone].nlist.keys()).difference(included_set)
				for node in lost_nodes:
					result_structure[zone].nlist.pop(node)
				for node in result_structure[zone].nlist.values():
					if node.num in included_set:
						node.norm()
		
		
		# Экспорт в файл или вернуть как структуру в случае если outname is None
		if hcnname is not None:
			export_to_mat_zone_bin_file(hcnname, result_structure)
			if verbose:
				print("\nFile {} has been written".format(hcnname))
			return None
		else:
			return result_structure


# Функция записывает файл в котором есть режимы с последовательностью моментов вермени
def create_regtime_file(file, sequences):
	if all(map(lambda a:isinstance(a, collections.abc.Iterable), sequences)):
		with open(file, 'w') as f:
			for n, i in enumerate(sequences, 1):
				f.write('r{n}_1 "r\t\t" {s}\n'.format(n=n, s=str(list(i))[1:-1].replace(',', '')))
	

# Метод определения абсолютного сжатия в узле
class CheckAbsoluteCompressMethod(enum.Enum):	
	ALL=0	
	NORMAL=1

# Функция извлекает необходимый тензор из фала типа hcn(hdf) в файлы типа tmp
# filelist - список имен(или путей) hcn(hdf) файлов
# nodefile - файл с номерами узлов которые необходимо извлеч(первая строка в этом файле комментарий). Елси значение None - извлекает все узлы.
# outdir - папка в которую извлекаются узлы
# verbose - выводить в консоль больше информации
# zeros - флаг который регулирует добавление нулевого момента времени к файлу tmp в первую строку
# tail - флаг который регулирует добавление нулевого момента времени к файлу tmp в последнюю строку
# listfile - флаг который создает в каждой папке материалов файл NODE.TXT со списком узлов в данной папке
# regtime - флаг регулирующий содание файла режимов (файл n.his в папке запуска скрипта)
# check_acompress - флаг который включает проверку узла на абсолютное сжатие
# check_acompress_method - метод определения полного сжатия в узле
# check_acompress_tolerance - предел значения при котором узел считается узлом с абсолютным сжатием
# compression_macro - флаг который регулирует запись макроса с узлами с абсолютным сжатием
# dont_extract - этот флаг регулирует если у нас нет необходимости записывать извлекаемые данные на диск а нужно к примеру извлечь только compression_macro
def extract_tensor_from_hcn(filelist, nodefile=None, outdir=None, verbose=False, zeros=False, tail=False, listfile=False, regtime=False, check_acompress=False, check_acompress_method:CheckAbsoluteCompressMethod=CheckAbsoluteCompressMethod.ALL, check_acompress_tolerance=0.0, compression_macro=False, dont_extract=False):
	global log
	if verbose:
		print('The actual list of files with data: {}'.format(str(filelist)[1:-1].replace("'",'')))
	if check_acompress_method == CheckAbsoluteCompressMethod.ALL:
		check_acompress_value = 7
	elif check_acompress_method == CheckAbsoluteCompressMethod.NORMAL:
		check_acompress_value = 4
	
	# Словарь из элемнтов:ключ - номер материала, значение - множество узлов в этом материале
	mat_node_set = {}

	if nodefile is None:
		nodeset = None
	else:
		nodeset = extract_nodes(nodefile)
		if verbose:
			print("Nodefile {}".format(nodefile))
			print("Length of nodefile:{}".format(len(nodeset)))
	if outdir is None:
		outdir = os.curdir
	curnum = 1
	
	temp = None 	# Шаблон вставки
	
	if regtime:
		regtimelist = []
	
	for file in filelist:
		with h5py.File(file, mode='r') as f:
			for mat in f.keys(): 
				# Определяем множество mat_node_set как ПЕРЕСЕЧЕНИЕ множества узлов в множестве nodeset и множества узлов в файле
				if nodeset is not None:
					mat_node_set[mat] = nodeset.intersection(set(map(int, f[mat].keys())))
				# Если nodeset отсутвсвует просто берем набор узлов в файле как множество mat_node_set
				else:
					mat_node_set[mat] = set(map(int, f[mat].keys()))
	
	if check_acompress:
		msz = {} # Словарь с зонами в которых будут множества словари с узлами в которых есть растягивающие напряжения больше критерия check_acompress_tolerance
		
		# Проверка на абсолютное сжатие
		for file in filelist:
			with h5py.File(file, mode='r') as f:
			# Для каждого материала в файле
				for mat in f.keys():
					try:
						cmsz = msz[mat]
					except KeyError:
						msz[mat] = set()
						cmsz = msz[mat]				
					for node in mat_node_set[mat]:
						if node not in cmsz:
							if numpy.max(numpy.max(f[mat][str(node)], 0)[1:check_acompress_value]) > check_acompress_tolerance:
								cmsz.add(node)
		for mat in msz.keys():
			if compression_macro:
				compression_dict = {i:0.0 for i in set.difference(mat_node_set[mat], msz[mat])}
				if len(compression_dict):
					if verbose:
						print("Write list of comperssion nodes in mat_compression_macro_{}.mac. Number of nodes in absolute compression : {}.".format(mat, len(compression_dict)))
					unite.generate_visulise_macro(compression_dict, "mat_compression_macro_{}.mac".format(mat), no_comment=True)
			if not dont_extract:
				mat_node_set[mat].intersection_update(msz[mat])
		
			
	# Для всех файлов в списке файлов hcn
	if not dont_extract:
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
					if curnum == 1:
						if len(mat_node_set[mat]):
							try:
								os.mkdir('{o}{sep}MATERIAL_{m}'.format(m=mat, sep=os.path.sep, o=outdir))
							except FileExistsError:
								pass
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
						if file == filelist[-1] and tail:
							of.write(FOUR_TEMPLATE.format(a = (20.0, 0.0, 0.0, 0.0, 0.0), num=n + 1, rname="NULL"))
						of.close()
			if regtime:
				regtimelist.append(range(curnum, n + 1))
			# Принимаем curnum n+1 для того чтобы момент времени из следующего файла пошел под сквозной нумерацией 
			curnum = n + 1
		if listfile:
			for mat in mat_node_set.keys():
				if len(mat_node_set[mat]):
					with open("{o}{sep}MATERIAL_{m}{sep}ZONA1.txt".format(m=mat, sep=os.path.sep, o=outdir), mode='w') as zf:
						zf.write("%i//\n" % len(mat_node_set[mat]))
						for i in sorted(mat_node_set[mat]):
							zf.write("%i\n" % i)
		if regtime:
			create_regtime_file('n.his', regtimelist)

# Класс RTable значительно упрощает чтение отдельных таблиц(записей) в файле типа rst
class RTable:
	TYP = {
	'i': 1,
	'q': 2,
	'd': 2,
	'f': 1
	}
	
	def __init__(self, file, name=None, size=None, typ=None):
		self.curpos = file.tell()
		if name is not None:
			self._name = name
		else:
			self._name = ''

		self._head = file.read(8)
		if size is None:
			size = struct.unpack('2i', self._head)[0]
		
		if typ is None:
			typ = struct.unpack('2i', self._head)[1]
			if typ == -2147483648:	# binary 1000 0000 0000 0000 0000 0000 0000 0000
				self._type = 'i'
			elif typ == 0:			# binary 0000 0000 0000 0000 0000 0000 0000 0000
				self._type = 'd'
			elif typ == 1073741824:	# binary 0100 0000 0000 0000 0000 0000 0000 0000
				self._type = 'f'
			else:
				raise TypeError('Unknown incorporated type')
		else:
			self._type = typ
		
		self._body = list(struct.unpack('{c}{t}'.format(c=int(size // RTable.TYP[self._type]), t=self._type), file.read(size * 4)))
		self._tail = file.read(4)
		if size != struct.unpack('i',self._tail)[0]:
			file.seek(self.curpos, 1)
			raise IndexError("{}-{}".format(size, struct.unpack('i', self._tail)[0]))
	
	@property
	def name(self):
		return self._name
	
	def __getitem__(self, item):
		return self._body[item]
	
	def __setitem__(self, item, value):
		self._body[item] = value
	
	def __iter__(self):
		return iter(self._body)
	
	def __len__(self):
		return len(self._body)
	
	def __bytes__(self):
		return self._head + struct.pack('{c}{t}'.format(c=len(self), t=self._type), *self._body) + self._tail

	def __str__(self):
		return str(list(enumerate(self._body)))[1:-1]
		
	def __repr__(self):
		return self._name + ":{}".format(str(self))

# Функция вставляющая необходимые значения, находящиеся в словаре damaged_dict rst
# rstfilename - файл rst который будет использован в качестве прототипа (для того чтобы взять координаты номеров узлов и прочее)
# new_rstfilename - вновь создаваемый файл
# damaged_dict - словарь с особой структурой в которой содержаться значения для записи (фактически поврежденность)
def get_damage_rst(rstfilename, new_rstfilename, damaged_dict, verbose=False):
	global TABLE_OF_ETYPES
	
	enod_table = {} # table of elements
	mat_node_table = {} #table of materials with nodes
	
	with open(rstfilename, 'rb', buffering=0) as f:
		# Тут происходит последовательная считка заголовков
		# согласно ANSYS Help help/ans_prog/Hlp_P_INT1_2.html (input that addr into the "Go to page" interactive form)
		info_header = RTable(f)
		rst_header = RTable(f)
		dof_header = RTable(f)
		node_table = RTable(f)
		elm_table = RTable(f)
		
		glbnnod = rst_header[48]
		glbnod_ptr = rst_header[49]
		gnod_table = None
		if glbnnod != len(node_table) and glbnod_ptr == 0:
			gnod_table = RTable(f)
		
		dsi_table = RTable(f)
		nsets = rst_header[8]
		start_point = dsi_table[0]
		if nsets == 1:
			end_point = info_header[26]
		else:
			end_point = dsi_table[1]
		
		time_table = RTable(f, typ='d')
		lsp_table = RTable(f)
		
		cyc_ptr = rst_header[16]
		cyc_table = None
		if cyc_ptr:
			cyc_table = RTable(f)
		
		ntran = rst_header[28]
		ntran_table = None
		if ntran:
			ntran_table = RTable(f)
		
		geo_header = RTable(f)
		ety_header = RTable(f)
		
		# Material configurations
		matnum = geo_header[13]
		for i in range(1, matnum + 1):
			mat_node_table[i] = MatZone(i)
			mat_node_table[i].nlist = {}
		# Read damage dict
		for num, item in damaged_dict.items():
			try:
				mat_node_table[num].nlist = item
			except KeyError:
				print('Material number {} hasnt met at the rst file'.format(num))
		
		ety_tables = []
		for i in ety_header:
			last_ety = RTable(f)
			ety_tables.append(last_ety)
			TABLE_OF_ETYPES[last_ety[0]] = (last_ety[1], last_ety[93])
		del last_ety
		
		rl_header = RTable(f)
		rl_tables = []
		for i in rl_header:
			if i!=0:
				rl_tables.append(RTable(f, typ='d'))
		
		maxcsy = geo_header[5]
		csy_header = None
		csy_tables = []
		if maxcsy:
			csy_header = RTable(f)
			for i in csy_header:
				csy_tables.append(RTable(f, typ='d'))
		
		loc_tables = []
		for i in node_table:
			loc_tables.append(RTable(f, typ='d'))
		
		eid_table = RTable(f, typ='q')
		
		# Configure elements
		element_tables = []
		for enum, i in zip(elm_table, eid_table):
			last_element = RTable(f)
			element_tables.append(last_element)
			try:
				_e = Element(enum, last_element[0], last_element[1], last_element[10:(10 + TABLE_OF_ETYPES[last_element[1]][1])])
			except IndexError:
				print(last_element)
			if _e is not None:
				enod_table[enum] = _e
		del _e
		del last_element
			
		another_info = f.read(start_point * 4 - f.tell())
		
		solution_header = RTable(f)
		solu_unused_info = []

		
		ptr_esl = solution_header[118] + solution_header[119]
		solu_unused_info.append(f.read((ptr_esl - len(solution_header) - 3) * 4))
		
		esl_table = RTable(f, typ='q')

		elemresults = []

		elems_walk = list(filter(lambda a: a[1]!=0, zip(elm_table ,esl_table)))
		elems_walk.sort(key=lambda a:a[1])

		for enum, elem in elems_walk:
			if elem == 0:
				continue
			try:
				curmat = enod_table[enum].mat
				nlist = enod_table[enum].nlist
			except KeyError:
				curmat = -100
				nlist = []
			# Это условие было добавлено для версии 16.2
			# В этой версии после данных о результатах в элементе существует 9 непонятных чисел, неотраженных в спецификации
			# Костыль добавляет их в общие результаты
			cur_elem_pos = (elem + ptr_esl + start_point) * 4
			if cur_elem_pos != f.tell():
				unknown_string = cur_elem_pos - f.tell()
				elemresults.append(f.read(unknown_string))
			elemresult = RTable(f)
			elemresults.append(elemresult)
			for inum ,item in enumerate(elemresult):
				if item > 0:
					ei = RTable(f, typ='f')
					if (ei.curpos - elemresult.curpos) // 4 == elemresult[2] and curmat != -100:
						for node, i in zip(nlist, range(0, len(ei), 6)):
							try:
								ei[i] = mat_node_table[curmat].nlist.get(node, 0.0)
							except KeyError:
								ei[i] = 0.0
							except AttributeError:
								raise
							ei[i + 1] = 0.0
							ei[i + 2] = 0.0
							ei[i + 3] = 0.0
							ei[i + 4] = 0.0
							ei[i + 5] = 0.0
					elemresults.append(ei)


		
		if f.tell() != end_point * 4: 
			print(elemresults[-20:])
			raise ValueError("{}-{}".format(f.tell(), end_point * 4))
		
		# CHANGE DSI TABLE
		for i in range(1, len(dsi_table) + 1):
			if dsi_table[i] == 0:
				break
			dsi_table[i] = 0
			
		
		# CHANGE NSETS
		rst_header[8] = 1
		
		# CHANGE TIME TABLE
		for i in range(1, len(time_table) + 1):
			if time_table[i] == 0.0:
				break
			time_table[i] = 0.0
			
		
		# CHANGE LSP TABLE
		for i in range(3,len(lsp_table) + 1):
			if lsp_table[i] == 0:
				break
			lsp_table[i] = 0
		
		# CHANGE INFO HEADER
		info_header[26] = end_point
		info_header[96] = end_point
		

						
		# Итоговая последовательность данных которая должна быть записана в rst файл
		result_sequence = (info_header, rst_header, dof_header, node_table, elm_table, gnod_table,
						dsi_table, time_table, lsp_table, cyc_table, ntran_table,
						geo_header, ety_header, *ety_tables, rl_header, *rl_tables,
						csy_header, *csy_tables, *loc_tables, eid_table, *element_tables,
						another_info, solution_header, solu_unused_info[0], esl_table)
		if verbose:
			print("Start the main write cycle...")
		with open(new_rstfilename, 'wb', buffering=0) as f:

						
			for item in result_sequence[:]:
				if item is None:
					continue
				elif isinstance(item, RTable):
					f.write(bytes(item))
				elif isinstance(item, bytes):
					f.write(item)
			
			for num, i in enumerate(elemresults):
				f.write(bytes(i))
				if verbose:
					if num % (len(elemresults) // 1000) == 0:
						print('Write results {}/{}\t\t'.format(num, len(elemresults)), end='\r')
			print('Write results {a}/{a}\t\t\n'.format(a=len(elemresults)), end='\r')


# Entry point function		
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
	# Вызов различных подпрограмм
	if args.extractor=='from-rst':
		base = None
		included_type = IncludedType.NON
		if args.elemfile is not None:
			included = extract_nodes(args.elemfile)
			included_type = IncludedType.ELEMENT
			if log:
				log.info('extract nodes which owned by elements which lie in {} file'.format(args.elemfile))
		elif args.nodefile is not None:
			included = extract_nodes(args.nodefile)
			included_type = IncludedType.NODE
			if log:
				log.info('extract nodes which contain in {} file'.format(args.nodefile))
		else:
			included = None
			if log:
				log.info('no exclusion for any item')
		if args.f:
			base = 6
		if args.a and log:
			log.info('border nodes between materials will be averaged by value')
		extract_tensor_from_rst(args.rstfile, args.hcnfile, base, included, included_type, args.a, args.v, args.sm, args.em)
	elif args.extractor=='from-hcn':
		extract_tensor_from_hcn(args.hcnfiles, args.nodefile, args.outdir, args.v, args.z, args.t, args.i, args.r, args.c or args.macro_compress, args.ca_method, args.ca_tol, args.macro_compress, args.macro_compress)
	elif args.extractor == 'to-rst':
		if unite:
			mdf = {}
			mnd = {}
			for material_num, dam_file in args.mdf:
				if mdf.get(material_num):
					mdf[material_num].append(dam_file)
				else:
					mdf[material_num] = [dam_file,]
			for mat, dam_files in mdf.items():
				# Использование модуля unite для извлечения знчение из нескольких макросов
				mnd[int(mat)] = unite.main(['-icf', '-jrd', '-s', '--com', '2', *dam_files]) 
				# При этом берутся максимальные значения если номера узлов содержаться в нескольких файлах
			get_damage_rst(args.rstfile, args.damfile, mnd, args.v)
		else:
			print(GLOBAL_UNITE_LACK_MESSAGE)
	else:
		print('print -h or --help key to show help message.')
	
if __name__=="__main__":
	main()