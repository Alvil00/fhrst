import sys
import argparse
import os

def arg_parse(arguments):
	parser = argparse.ArgumentParser(description = "unite visualise macro script for crp")
	parser.add_argument('--com',type=int,default=None,help="commented start lines number")
	parser.add_argument('--out',type=str,default=None,help="output filename")
	input_file_group = parser.add_mutually_exclusive_group()
	input_file_group.add_argument('-icf', action="store_true", help="input as dof conc macro file.Default option")
	input_file_group.add_argument('-ibf', action="store_true", help="input as bf macro file")
	input_file_group.add_argument('-idf', action="store_true", help="input as damage file")
	output_file_group = parser.add_mutually_exclusive_group()
	output_file_group.add_argument('-ocf', action="store_true", help="output as dof conc macro file. Default option")
	output_file_group.add_argument('-obf', action="store_true", help="output as bf macro file")
	output_file_group.add_argument('-odf', action="store_true", help="output as damage file")
	if __name__ != '__main__':
		output_file_group.add_argument('-jrd', action="store_true", help="dont output as file, just return united dict")
	parser.add_argument('-s',action="store_true", help="silence mode")
	parser.add_argument('-op', action="store_true", help="pipe mode")
	parser.add_argument('--method', default='max', choices=('max', 'min'), type=str, help='metod for compartion: chose max value or min')
	parser.add_argument('-no-output-comment' , action='store_true', help='no output file comment', dest='noc')
	parser.add_argument('files',metavar='macro_file',nargs='+',type=str,help="Macro files")
	args = parser.parse_args(arguments)
	if __name__ == '__main__':
		args.jrd = False
	
	if not(args.icf or args.ibf or args.idf):
		args.icf = True
		
	if not(args.ocf or args.obf or args.odf or args.jrd):
		args.ocf = True	
		
	return args

def read_nodes_from_node_damage_file(filename=None, format=(int,float), com_str=2 , s=False):
	if not filename:
		filename = "NODE.txt"
	# this list save node structure
	# each element of the list has format [node number,node zona] 
	# and maybe some another values included in format argument of function
	node_dict = dict()
	c = list()
	# void list to apply format logic to standalone values in line
	with open(filename, "r") as node_file:
		for num, i in enumerate(node_file):
			# first two strings in standard NODE.txt file which has created by ansalt program
			# consist from comments
			# then we need to skip them
			if num >= com_str:
				buf_str_list = i.split() # we split the node_str to standalone values 
				if len(buf_str_list)<len(format):
					for q in range(len(format)-len(buf_str_list)):
						buf_str_list.append("0")
				for format_num,format_cast in enumerate(format): # and apply to all of them appropriate cast
					c.append(format_cast(buf_str_list[format_num]))
				node_dict[c[0]] = float(c[2]) # append number of node and node group (and maybe some another values included in format argument) in a list:new version
				c = list() # clear list which we appended to returnable list
	return node_dict
	
def read_nodes_from_visualse_macro(filename = None, vistype=None, com_str = 4, s = False):
	if filename is None:
		filename = 'VIS.mac'
	if vistype == None:
		vistype = "DOF"
	node_dict = {}
	with open(filename,'r') as f:
		for num,i in enumerate(f):
			if num>=com_str:
				if i == "/go" or not i.strip():
					break
				node = i.split(",")
				if vistype=="DOF":
					node_dict[int(node[1])]=float(node[4])
				elif vistype=="BF":
					node_dict[int(node[1])]=float(node[3])
				else:
					if not s:
						print("unknown format")
					raise Error
	return node_dict

DOF_TEMPLATE = "DNSOL, {node_num},CONC,, {node_damage}\n"
BF_TEMPLATE =  "BF, {node_num},TEMP, {node_damage}\n"

def generate_visulise_macro(node_dict, filename= None, topipe = False, vistype = None, no_comment=False):
	global DOF_TEMPLATE
	global BF_TEMPLATE
	
	if filename is None:
		filename = 'VIS.mac'
	if vistype == None:
		vistype = "DOF"
		
	if vistype == "DOF":
		v_template=DOF_TEMPLATE
	elif vistype == "BF":
		v_template=BF_TEMPLATE
		
	with open(filename,"w") as mac_file:
		if not no_comment:
			mac_file.writelines("/nopr\n")
			mac_file.writelines("/com,--- This macros has been created by Unite Visualise Macro routine ---\n")
		if vistype=="DOF":
			mac_file.writelines("/GRAPHICS, FULL\n")
			mac_file.writelines("DOF,CONC\n")
		for i in node_dict.keys():
			mac_file.writelines(v_template.format(node_num=i,node_damage=node_dict.get(i)))
		if not no_comment:
			mac_file.writelines("/go")

def generate_node_damage_file(node_dict, filename=None, topipe = False):
	if filename is None:
		filename = 'NODE_DAMAGE.TXT'
	node_damage_list = list(zip(node_dict.keys(),node_dict.values()))
	node_damage_list.sort(key = lambda a:a[1],reverse = True)
	template = '{:15}\t0\t{:15}\n'
	if not topipe:
		with open(filename,'w') as f:
			for i in node_damage_list:
				f.write(template.format(i[0],i[1]))
	else:
		for i in node_damage_list:
			print(template.format(i[0],i[1]))
			
def main(arguments=None):

	if arguments is None:
		arguments = sys.argv[1:]

	args = arg_parse(arguments)
	
	if args.idf:
		ifunc = read_nodes_from_node_damage_file
		if args.com is None:
			args.com = 0
		ai = ((int,int,float),args.com,args.s)
	elif args.ibf:
		ifunc = read_nodes_from_visualse_macro
		if args.com is None:
			args.com = 2
		ai = ("BF",args.com,args.s)
	elif args.icf:
		ifunc = read_nodes_from_visualse_macro
		if args.com is None:
			args.com = 4
		ai = ("DOF",args.com,args.s)
	
	if args.method == 'max':
		cfunc = lambda a, b: a > b
	elif args.method == 'min':
		cfunc = lambda a, b: a < b

	node_list = []
	
	for file in args.files:
		try:
			node_list.append(ifunc(file,*ai))
		except FileNotFoundError as fe:
			if not args.s:
				print("{}: {}".format(fe.args[1],fe.filename))
			quit()
			
	united = {}
	for num,i in enumerate(node_list):
		if num>0:
			intersection_key_list = list(set(united.keys()).intersection(set(i.keys())))
			for j in intersection_key_list:
				if cfunc(united.get(j), i.get(j)):
					i[j]=united[j]
		united.update(i)
	if args.odf:
		ofunc = generate_node_damage_file
		ao = (united,args.out,args.op)
	elif args.obf:
		ofunc = generate_visulise_macro
		ao = (united,args.out,args.op,"BF",args.noc)
	elif args.ocf:
		ofunc = generate_visulise_macro
		ao = (united,args.out,args.op,"DOF",args.noc)
	elif args.jrd:
		return united
		
	ofunc(*ao)
	
if __name__=="__main__":
	main()