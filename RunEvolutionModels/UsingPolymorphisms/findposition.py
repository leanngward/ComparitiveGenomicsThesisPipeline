import sys


filename = sys.argv[1]
infile = open(filename,'r')
line = infile.readline()
#list = ['7599','7600','7601','7602','7603','7604','7605','7606','7607','7608','7745','7751','7769','7770','7773','7841','7871','7874','8198','9748','9880','14454','14531','14550','14557','14588','14594','30597']
list = ['306107','306108','306109','306110','306111','306112','306113','306114']
count = 1
for items in line:
	if str(count) in list:
		print("count:")
		print(count)
		print("nuc:")
		print(items+"\n")
	count += 1

infile.close()
