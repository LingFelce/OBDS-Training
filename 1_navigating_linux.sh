#Use PUTTY or Terminal in X2Go

#Log into CGAT system
ssh cgatui.imm.ox.ac.uk

#Log into head node
ssh cgath1

#type command shows how command is interpreted
type <command> 

#which command displays which executable programme will be executed (and where it's located)
which <command>

#displays a command's manual page
man <command>

#displays command's info entry
info <command>

#displays brief description of command
whatis <command>

#print working directory (where you are)
pwd

#list directory contents with options
ls -h #human readable
ls -a #all files
ls -F #will append / if directory
ls -s #size
ls -l #long format

#change directory
cd

#print a brief escription of file's content
file <filename>

#print file contents to standard output
cat <filename>

#print first 20 lines
head -20 <filename>

#print last 20 lines
tail -20 <filename>

#display output in terminal one page at a time
more <filename>

#like more but allows backward movement. Also works for zipped files
less <filename>

#opens nano terminal text editor
nano <filename>

#counting characters, lines or words
wc -m <filename> #character counts
wc -l <filename> #print new line counts
wc -w <filename> #print word count
wc --help #display help message

#create new empty file
touch <filename>

#create new directory (can do multiple ones)
mkdir <directory>

#copy files
cp <file1> <file2> #copy file to new file
cp <file1> dir/ #copy file to new directory

#move files
mv <file1> <file2> #rename file
mv <file1> dir1/ #move file to new directory

#rename multiple files using pattern matching
rename 's/old/new/' <files>
rename 's/perl/pl/' *.perl

#removing files
rm file1 #remove file
rm *.perl #remove multiple files with same file extension
rm -r dir1 #remove directory and contents
rmdir dir1 #remove empty directory

#searching for files
locate <filename>
locate text.png
find ~ -type f -name ".JPG" #search in ~ home . current / root (everywhere) search for .JPG files
find . -name "*.tsv" -exec wv -l {} #execute command word count on each file found, brackets signify where to put current file name

#comparing files
comm file1.txt file2.txt #compares files line by line. 1st col lines unique to file1, 2nd col lines unique to file2, 3rd col shared lines
diff file1.txt file2.txt #more advanced tool

#file compression
gzip <file> #compress file in current directory
gunzip <file> #unzip
gunzip -c <file> #decompress to standard output
zcat <file> #print compressed file to stdout
zless <file> #less compressed file

#hard and symoblic links
ln file1 name2 #create hard link
ln -s file1 link1 #create symbolic link

#file transfer
wget http://url #used to mirror websites

#.bashrc is shell script run every time user opens new shell. Set environment variables, aliases, virtual environments
#.inputrc customise how keystrokes intpreted in terminal

#bash shortcuts
# CTRL-c Abort current command
# CTRL-z Pause current foreground process
# CTRL-l Clear the screen
# CTRL-a Go to start of line
# CTRL-e Go to end of line
# CTRL-u Cut from start of line
# CTRL-k Cut to end of line
# CTRL-r Search history
# CTRL-d Logout (also exit)
# Up arrow – access previous commands
# Tab – autocomplete (will prompt if ambiguous)

#changing file permissions
chmod ug+rwx file #user and group have read, write and execute permissions
ls -l file #check permissions for file

chmod 770 file #same as above but using octal notation instead of alphabetical

d	r	w	x	r	w	x	r	w	x
Owner	Group	Other
Directory	Read	Write	Execute	Read	Write	Execute	Read	Write	Execute

drwxr-xr-x
A folder which has read, write and execute permissions for the owner, but only read and execute permissions for the group and for other users.
-rw-rw-rw-
A file that can be read and written by anyone, but not executed at all.
-rw-r--r--
A file that can be read and written by the user, but only read by the group and everyone else.

HOME / BLOG / RUNNING A WEBSITE / UNIX AND LINUX / UNDERSTANDING FILE PERMISSIONS
Understanding File Permissions
24 APRIL 2001 / 14 COMMENTS

One of the hardest things for the beginner webmaster to get to grips with is how to use chmod correctly to set permissions on files on UNIX and Linux web servers. You need to set the correct permissions on CGI scripts when you install them, to stop those dreaded “500 Server Error” messages.

In this tutorial, I’m going to explain the concept of permissions, and show you how to set permissions using your FTP program or via SSH.


What are permissions?
On a UNIX web server, every single file and folder stored on the hard drive has a set of permissions associated with it, which says who is allowed to do what with the file. Every file (and folder) also has an “owner” and a “group” associated with it. If you created the file, then you are usually the owner of that file, and your group, or the group associated with the folder you created the file in, will usually be associated with that file.

Who can do stuff?
There are three types of people that can do stuff to files – the Owner of the file, anyone in the Group that the file belongs to, and Others (everyone else). In UNIX, these 3 types of people are referred to using the letters U (for Owner, or User in UNIX-speak!), G (for Group), and O (for Others).

What stuff can you do?
There are three basic things that can be done to files or folders:

You can read the file. For folders, this means listing the contents of the folder.
You can write to (change) the file. For folders, this means creating and deleting files in the folder.
You can execute (run) the file, if it’s a program or script. For folders, this means accessing files in the folder.
What do all these funny letters and numbers mean?!
That’s the basics of permissions covered. As you can see, there’s not much to them really!

The confusion often occurs when you have to start actually setting permissions on your file server. CGI scripts will tell you to do things like “chmod 755” or “Check that the file is executable”. Also, when you use FTP or SSH, you’ll see lots of funny letters next to the files (such as rwxrw-rw-). We’ll now explain what all these hieroglyphics mean!

When you FTP to your web server, you’ll probably see something like this next to every file and folder:

Attributes list
This string of letters, drwxrwxrwx, represents the permissions that are set for this folder. (Note that these are often called attributes by FTP programs.) Let’s explain what each of these letters means:

# d	r	w	x	r	w	x	r	w	x
# Owner	Group	Other
# Directory	Read	Write	Execute	Read	Write	Execute	Read	Write	Execute

# drwxr-xr-x
# A folder which has read, write and execute permissions for the owner, but only read and execute permissions for the group and for other users.
# -rw-rw-rw-
# A file that can be read and written by anyone, but not executed at all.
# -rw-r--r--
# A file that can be read and written by the user, but only read by the group and everyone else.

chmod 777 # is the same as rwxrwxrwx
chmod 755 # is the same as rwxr-xr-x
chmod 666 # is the same as rw-rw-rw-
chmod 744 # is the same as rwxr--r--

#managing Linux processes
bg #run paused process in background (or <command> & to start process running automatically in background)
fg #bring background process to foreground
ps #status of processes for user
kill <PID> #kills process with specific process ID
top #details on active processes
df #free hard disk space
du #directory usage space
free #free RAM


