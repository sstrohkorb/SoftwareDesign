dna = 'ACTGAGTCATTTTTTAAG'

for i in range(0, len(dna), 3):
  print dna[i:i+3]

# output the index where the desired number is located in the list
def get_index(list1, desired, list2):
  for i in range(len(list1)):
    #print list1[i]
    if list1[i] == desired:
      #print "Happiness"
      return list2[i]


list1 = [1,2,3]
list2 = ['a','b','c']
desired = 3
print get_index(list1, desired, list2)


input_a = 2
print "Input = " + str(input_a)
expected_output_a = 'b'
print "Expected Output = " + str(expected_output_a)
actual_output_a = get_index(list1, input_a, list2)
print "Actual Output = " + str(actual_output_a)

