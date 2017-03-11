import sys
import random

final = []
numbers = random.sample(range(0, sys.maxsize), 20)
for num in numbers:
	final.append('{0:078b}'.format(num))

file = open('msg.txt', 'w+')
for binary in final:
	file.write(binary + '\n')
file.close()
