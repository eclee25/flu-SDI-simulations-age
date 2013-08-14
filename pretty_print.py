# module to print data structures nicely

def print_list(l):
	for value in l:
		print "%s" % value

def print_dictionary(dic):
	for key, value in dic.items():
		print "%s,%s" % (key,value)

def print_list_to_file(l, filename):
	fwriter = open(filename, 'w+')
	for value in l:
		fwriter.write("%s\n" % value)
	fwriter.close()

def print_dictionary_to_file(dic, filename):
	fwriter = open(filename, 'w+')
	for key, value in dic.items():
		fwriter.write("%s,%s\n" % (key,value))
	fwriter.close()

def print_dict_to_file_2k3v(dic, filename):
	fwriter = open(filename, 'w+')
	for k, v in dic.items():
		fwriter.write('%s,%s,%s,%s,%s\n' % (k[0], k[1], v[0], v[1], v[2]))
	fwriter.close()

def print_dictlist_to_file(dic, filename):
	fwriter = open(filename, 'w+')
	for k, v in dic.items():
		str_v = str(v)[1:-1]
		fwriter.write('%s\n' % (str_v))

def print_OR_to_file(dic, filename):
	fwriter = open(filename, 'w+')
	for k, v in dic.items():
		str_k = str(k)
		str_v = str(v)[1:-1]
		fwriter.write('%s,%s\n' % (str_k, str_v))