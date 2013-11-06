# module to print data structures nicely

# dependent packages
import zipfile

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
		str_v.replace(" ", "") # remove all internal whitespaces to save disk space
		fwriter.write('%s\n' % (str_v))
	fwriter.close()

def print_OR_to_file(dic, filename):
	fwriter = open(filename, 'w+')
	for k, v in dic.items():
		str_k = str(k)
		str_v = str(v)[1:-1]
		fwriter.write('%s,%s\n' % (str_k, str_v))
	fwriter.close()

def print_OR_time_to_file(dic, filename, beta):
	fwriter = open(filename, 'w+')
	for k, v in dic.items():
		if k[0] == beta:
			str_k = str(k)
			str_v = str(v)[1:-1]
			fwriter.write('%s\n' % (str_v))
	fwriter.close()

def print_sorteddlist_to_file(dic, filename, numsims):
	fwriter = open(filename, 'w+')
	for k in xrange(numsims):
		str_v = str(dic[k])[1:-1]
		str_v.replace(" ", "") # remove all internal whitespaces to save disk space
		fwriter.write('%s\n' % (str_v))
	fwriter.close()

def compress_to_ziparchive(zipname, filename):
	with zipfile.ZipFile(zipname, 'a') as zf:
		zf.write(filename, compress_type = zipfile.ZIP_DEFLATED)

	
	
	
			
			