# Fibonacci series:
# the sum of two elements defines the next
a, b = 0, 1

# read n from stdin
n = int( input( "Enter maximum:" ))
print( "maximum is", n )
while a < n:
	print(b)
	a, b = b, a+b


a, b = 0, 1
count = 0;
while count < n:
	count += 1
	print( "count,value:", count, b )
	a, b = b, a+b
