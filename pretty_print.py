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

def print_dictionary_to_file(dic, filename):
	fwriter = open(filename, 'w+')
	for key, value in dic.items():
		fwriter.write("%s,%s\n" % (key,value))

