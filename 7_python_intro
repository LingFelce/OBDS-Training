#7_Introduction to Python 23rd September 2019

#load X2Go. In terminal, load environment first then type python to start python console - should be version 3.6.7
#have to be in terminal + environment to load Spyder = integrated development environment for writing, testing, debugging, code analysis 
#Within Spyder, use F9 to run current line in text editor in console. Use control + enter to run code cells (like R chunks, need to set up separately)
#Jupyter Notebook - like .Rmd in R Studio, combine plots and code. JupyterLab new version. Runs in browser, shareable.

#github repository /ifs/obds-training/lingf/obds/devel/OBDS_Training_Sep_2019

#variables in Python must start with letter or underscore, case sensitive (usually use lower case)
#numbers (integer, float); strings (can use single, double or triple quotes), Boolean.
#operators (assignment, arithmetic, comparison, logical) == equal to, != not equal to
#if and loop (for and while)
#use if first, then if not true go on to elif (if use just if then will run all ifs)

#flow control - conditional, logical test
number = 23
guess = int(input('Enter an integer : '))
if guess == number:
    # New block starts here
    print('Congratulations, you guessed it.')
    print('(but you do not win any prizes!)')
# New block ends here
elif guess < number:
    # Another block
    print('No, it is higher than that')
    # You can do whatever you want in a block ...
else:
    print('No, it is lower than that')
    # you must have guessed > number to reach here

#flow control - iteration, loop

for i in range(1, 5):
print(i)
running = True
while running:
    guess = int(input('Enter an integer : '))
    if guess == number:
        print('Congratulations!')
        # this causes the while loop to stop
        running = False
elif guess < number:
    print('No, it is higher than that.')
else:
    print('No, it is lower than that.')

#python written in code blocks (tab/four spaces indentation), not brackets like R. Style guide - PEP8

#use def keyword to give name to block of statements (functions - like having separate functions.R file)
def say_hello():
    print('hello world')

say_hello() #call function

#return data from function (then can print)

#assign global variables in outermost code block, so accessible outside functions

#data structures (list (denoted using []), tuple (denoted using (), unchangeable list - no append, no sort), dictionary (denoted by {}, key-value pairs) - are all sequences)
#probably won't use tuple that often

len(listname) #find out number of items
listname.append('item') #add item to list
listname.sort() #sort list
del listname[0] #delete first item from list #[] used for slicing 

#slicing lists
my_list=['p','r','o','g','r','a','m']
print(my_list[2:5]) #elements 3rd to 5th
print(my_list[:-3]) #elements beginningto 4th
print(my_list[:2]) #elements 1st to 3rd
print(my_list[5:]) #elements to end
print(my_list[:]) #elements beginning to end
#output 
['o', 'g', 'r']
['p', 'r', 'o', 'g']
['p', 'r']
['a', 'm']
['p', 'r', 'o', 'g', 'r', 'a', 'm']

#nest lists - data structures within data structures to create high dimensional structures
nested=["hello", 2.0, 5, [10,20]]
nested[3][1]
#output 20

#by default don't copy lists as can be very large!
mylist=shoplist #mylist is another name for shoplist, any changes made will be the same
mylist=shoplist[:] #make a copy by doing a full slice

#dictionaries like an address-book, keys associated with values, keys must be unique and immutable, not ordered (randomised for security reasons)

ab = {
    'David': 'david.sims@imm.ox.ac.uk',
    'Charlie': 'charlotte.george@imm.ox.ac.uk',
    'Spammer': 'spammer@hotmail.com'
}
print("David’s email is", ab['David'])
print('There are ‘, len(ab), ‘ contacts’)
# Returning a value
email = ab.get(‘David’, 0)
# Returns 0 if not found
# Deleting a key-value pair
del ab['Spammer']
# Adding a key-value pair
ab['Guido'] = 'guido@python.org'
# Searching the keys
if 'Guido' in ab:
    print("Guido's address is", ab['Guido'])

#make shallow copy (does not copy nested data structures)
ab2=dict(ab)
ab2=ab.copy()
#deep copy (copies nested data structures)
ab2=copy.deepcopy(ab)

#read and write files
#open for 'w'riting
f=open('name.txt', 'w')
f.write(sequence) #write text to file
f.close() #close file

#open for 'r'eading
with open('name.txt','r') as f: #use with for opening file, automatically closes file when get to end of code block, good if code crashes
    for line in f:
        print(line)

#string manipulation - a string is an object and a sequence, are sliceable

#code organisation levels - functions > scripts > modules > packages (pandas) > libraries (python standard library)

#python scripts - designed to be run on command line so needs a main function at top e.g #!/usr/bin/env python (where to start running code). Script layout:
"""documentation docstring"""
#import useful existing code
#global variable definitions
#functions
#main function

#modules - set of functions with common theme, modules meant to be imported, no main function
#package - folder with set of files with software, documentation, tests (scripts and modules)
#library - collection of packages, exposes Application Programming Interface (API), can use to develop new software on top of it
#python standard library docs.python.org/3/library/index.html

#ALGORITHMS - input, algorithm (instructions), output. Can use flowchats to design, pseudocode (steps written in plain English)

#finding max number in a list - create a loop. 1st number = maxvalue. If new number is bigger than maxvalue, then that becomes new maxvalue.
#keeps going and comparing all numbers on list
numbers=['26','54','93','17','77','31','44','55','20']
maxvalue=numbers[0]
for i in numbers: #i is number in list
       if maxvalue<i:
           maxvalue=i
print(maxvalue)

#selection sort
numbers=['26','54','93','17','77','31','44','55','20']

#write a function to find the smallest value in list, return number and position
def min_value_func(num_list):
    min_value=num_list[0] #set min_value to be 1st number in list
    for i in num_list:
        if min_value > i: #if number in list is smaller than min_value
            min_value = i #then number becomes new min_value
    return min_value, num_list.index(min_value) #print min_value and positions

print(min_value_func(numbers))

for i in numbers:
    idx = numbers.index(i) #idx is position of number i
    curr_min = min_value_func(numbers[idx:]) #curr_min tuple containing min_value and position
    numbers[idx] = curr_min[0] #move min_value into lowest currently available position on list
    numbers[curr_min[1]+idx] = i #shuffle numbers in list rather than overwriting
    print(numbers)

#bubble sort
numbers =['26','54','93','17','77','31','44','55','20']

n = len(numbers)

for j in range(0,n-2): #outer loop
    offset = len(numbers)-j #which subset of outer loop you want to run inner loop to not repeat
    for i in range(0, offset-1): #finds all positions in list -1 (otherwise get error list index out of range)
        if numbers[i] > numbers[i+1]: #if number is larger than number in position above
            temp_var = numbers[i] #make temporary variable
            numbers[i] = numbers[i+1] #swap numbers around but don't overwrite
            numbers[i+1] = temp_var
    print(numbers)

#find GC content
sequence = "ATGCATGCATGC"
GC_content = sequence.count("G") + sequence.count("C")
GC_percentage = (GC_content / len(sequence))*100
print(GC_percentage, "%",sep="") #no separation between number and %

def gc_count(sequence):
    GC_content = sequence.count("G") + sequence.count("C")
    GC_percentage = (GC_content / len(sequence))*100
 
    print(GC_percentage, "%",sep="")

with open('/ifs/obds-training/lingf/obds/devel/OBDS_Training_Sep_2019/Homo_sapiens_CEBPB_sequence.fa.txt', 'r') as fasta:
    header = "" #null
    sequence="" #null
    first_line = True
    for line in fasta:
        if line.startswith (">"): #find line that starts with >
            if first_line:
                first_line = False #found first line, no more first lines after
            else: 
                print(header)
                gc_count(sequence)
            header = line.strip() #call it a header, removes white space
            sequence="" #resets sequence so can recognise next block
          
        else:
            sequence += line #append so that looks at all lines
    print(header) #last chunk
    gc_count(sequence) #last chunk


#DEBUGGING - syntax errors, semantic errors, algorithmic errors, complex errors
#avoid errors - write legible code, (write tests), write incrementally, reuse code, write documentation, use a debugger (pdb - integrated into spyder, blue symbols at top)
#extensive logging (use logging facility - python logging module part of standard library) - output what you are going to do, output what you've done, output anything unexpected and of interest
#debugging - to use pdb in spyder, can set breakpoint (double click margin so red dot appears) so runs debugger until that point
