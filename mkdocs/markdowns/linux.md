# Shell and Linux

> Author: Hanrui Bai, 2021-2-15

## About Linux
- The complete Linux operating system includes the kernel and outer layer applications. Different Linux distributions use the same kernel and different outer layer applications. 
- Common Linux distributions are [Ubuntu](), [red hat](), [SuSE](), [Gentoo](), [CentOS]() and [Debian](). 
- Shell is the interface between the linux kernel and the user.

## Get Linux

1. Install Ubuntu in a virtual environment.<Font color = red> Not suggested to novices.</Font>
2. Buy a Ubuntu.
3. Apply for an account at the Institute’s Supercomputer Center. <Font color = red> This is the most common method.</Font>

#### Log into Linux sever from Windows

1. Download a Secure Shell (ssh) terminal.
Common terminal:
    - Xshell: [https://www.xshellcn.com/](https://www.xshellcn.com/)
    - PuTTY: [http://www.putty.org/](http://www.putty.org/)

2. Install the shell terminal.

3. Create a new session.

4. Choose the ssh protocol and enter the host IP, port, user ID and password.

5. Click `connect`.

#### Log into Linux sever from Linux or Mac

Type command in terminal: `ssh [host IP:port]@[user ID]`. Then enter the password.

## Basic knowledge of Linux

1. The absolute path is identified by a forward slash (`/`), and the structure of directory is tree structure.

2. How to use linux command: enter the command at the shell prompt: `$`. Command like: `[Path]/command [-option parameter]  [file|directory]`.
For example: `~/miniconda3/envs/RNA/bin/fastp -i /home/sysbio/SRR8467686.fastq.gz`

3. Enter a single dot (`.`) to indicate the current directory, and enter a double dot (`..`) to indicate the parent directory of the current directory.
For example: `cd ..`

4. Use the man command to view the operation manual of all required commands. Command like: `man [command name]`.
For example: `man cp`

5. Filter: The filter refers to the name of the specified file and is used after many commands. Add a filter after the ls command, then only the information of the file will be displayed.
For example: `ls -l my_file`
This command specifies the relevant information of the output file my_file.
<Font color = red> When writing a filter, you can use a question mark (`?`) to represent a character, an asterisk (`*`) to represent zero or more characters, and use the tab key (`Tab`) to quickly complete the file name or directory first name.</Font>

6. The `which` command is used to find the path of the command. Command like: `which grep`
For example: `which fastp`

## About permission

#### View the file permission

When you use `ls -l` command, the permissions will be marked at the beginning of the file, like: drwxrwxrwx. The first letter represents the type of file: `d` means directory, `-` means file, and `l` means soft link. The remaining nine letters are grouped into three, representing the permissions for the owner, the owner`s` group, and other groups. `r` means read-only. `w` means allowing user to modify the content of the file. `x` means allowing user to execute the file as a program. `-` means having no such permission.

#### Modify permissions: chmod command

<Font color = red>You must have the authority to operate files</Font>
Command like: `chmod (user permissions) (group permissions) (other permissions) file`
Permission is expressed in octal code: r=4, w=2, x=1
For example: `chmod 755 test.txt`
The result is `-rwxr-xr-x test.txt`

## Command about directory

#### Check files which are in the directory list

1. Basic list function: ls command
  Command like: ls various parameters directory name
    - `-F` Distinguish between files and directories
    - `-R` Recursive option to list files in subdirectories contained in the current directory
    - `-l` Produces a long list of output results, including information about each file
    - `-d` Lists only the information of the directory itself, not its contents
    - `-i` The inode number of the file, which is the unique identifier of the file
    - `-a` Display both hidden files and ordinary files
    - `-t` Sort by time modified
    - `-h` Show file size easy to read for human
<Font color = red>The parameters can be combined and written. For example: `ls -ltrh`.</Font>

2. The way to output the tree list: tree tool
   The tool needs to be installed by yourself. Command like: `tree ./`.

#### Create a directory: mkdir command

Command like: `mkdir [directory name]`
To create multiple directories and subdirectories, you can use the `-p` parameter. For example: `mkdir -p New_file/work/file1`.
<Font color = red>The `-p` parameter in the mkdir command can create missing parent directories as needed.</Font>
   
#### Delete directory: rmdir command

Command like: `rmdir [-option parameter] [directory name]`
- `-No` Parameters after rmdir delete empty directories
- `-rf` Delete all contents in the directory
- `-i` Ask a question before deleting

<Font color = red>Because Linux does not have a recycle bin, you must add the -i parameter when deleting to confirm whether the deleted content is correct.</Font>

#### Switch the directory: cd command

Command like: `cd directory`
For example: `cd /usr/bin`

#### View the current absolute path of current directory

Command like: `pwd`
  
## Command about file

#### Create a file: touch command

Command like: `touch [-option parameter] [file name]`
- Create a new file without adding parameters after touch, if the file already exists, change the modification time.
- `-a` Change the access time of an existing file
  
#### Delete files: rm command

Command like: `rm -i [file name]`
<Font color = red>Because Linux does not have a recycle bin, you must add the `-i` parameter to avoid errors when deleting.</Font>
  
#### Copy files: cp command

Command like: `cp [-option parameter] [file name]`
- `-i` Force to ask if you need to overwrite existing files
- `-R` Copy the entire directory recursively
<Font color = red>To avoid errors, it is recommended to add the -i parameter</Font>

## Rename and move files: mv command

1. Rename
  Command like: `mv [original file name] [new file name]`
  
2. Move
  Command like: `mv [file name in the original path] [target new path]`
  
3. Rename and move operations can be performed at the same time
  For example: `mv /home/picture/book /home/file/cook`. Move the `book` file in the `picture` folder to the `file` folder, and rename it to `cook`.
  
## View files

#### View file type: file command

Command like: `file [file name]`
You can check the file type, character encoding method, whether the file can be run, etc.
 
#### View the entire file content

1. `cat` command: display all data in the text file
Command like: `cat [-option parameter] [file name]`
    - No parameters after cat means display content.
    - `-n` Display content after adding line number.
    - `-b` Displays only after adding line numbers to text content.
    - `-T` Does not display tabs, replace tabs with ^T to display.

2. `more` command: display the content of the text file, but it will stop after each page is displayed.

3. `less` command: display the content of text files and support more advanced interaction.

### View a part of file content

1. `head` command: view the beginning of the file
Command like: `head -n file name`
`n` is the number of rows displayed

2. `tail` command: view the end of the file
Command like: `tail -n file name`
`n` is the number of rows displayed

### Upload files, download files

1. `rz`: upload file
2. `sz [file name]`: download file

## Command about process

#### Probe the process

1. `ps` command
command like: `ps [-option parameter]`
Symbol description:
  - `PID`: Process ID.
  - `TTY`: Terminal device when process starts.
  - `TIME`: Cumulative CPU time required by the process.
  - `CMD`: The name of the command used.

2. top command
Command like: `top`
Symbol description:
  - `load average`: The three values are the average load of the last 1min, 5min, and 15min. The larger the value, the larger the load, and the more than 2 indicates the system is busy.
  - `PR`: Process priority
  - `NI`: Moderate value of process
  - `VIRT`: The total amount of virtual memory occupied by the process
  - `RES`: The total amount of physical memory occupied by process.
  - `MEM`: The ratio of memory used by the process to available memory.
  - `S`: Process state (`D`: interruptible sleep state; `R`: running; `S`: sleeping; `T`: tracking or stopping state; `Z`: rigid state).

3. The difference between ps and top command
The ps command displays information at a specific time point, and the top command displays real-time information.

## End the process

1. kill command
command like: `kill PID` or `kill -s [process signal]`

2. killall command
Command like: `killall [process name]`
<Font color = red>Use wildcards carefully in the killall command</Font>

## Command about task

#### Execute tasks to the background

1. `[Path]/command [-option parameter] [file|directory] &`
2. `nohup [path]/command [-option parameter] [file|directory] &`
3. screen command, command like :
  - Create a screen: `screen -dmS screen_test`
  - View the screen: `screen -list`
  - Connect to the screen: `screen -r screen_test`

#### Common task management commands:

1. `jobs`: View tasks, return task number n and process number
2. `bg %n`: Transfer task number n to the background
3. `fg %n`: Transfer task number n to the foreground
4. `ctrl+z`: Suspend the current task
5. `ctrl+c`: Stop the current task
6. `kill -n task`: End the task

## Command about disk space

#### Mount the device: mount command

Command like: mount parameter file device type device to be mounted target location
File device types are:
- `-vfat`: Windows long file system
- `-ntfs`: An advanced file system widely used in Windows
- `-iso9660`: Standard CD-ROM file system

#### Uninstall the device: umount command

Command like: `umount location/target device`

#### Explore the disk situation

1. df command
Command like: `df [-option parameter] [target disk]`

2. du command
Command like: `du [-option parameter] [target disk]`
   - `-c` Displays the total size of all listed files
   - `-h` Display in a user-readable way
   - `-s` Displays the total of each output parameter
   
3. The difference between df and du commands
The df command displays the disk usage, the du command displays the disk usage of each file

4. Disk partition: fdisk command
Command like: `fdisk -l [disk name]`

5. Disk formatting: mkfs command
Command like: `mksf -t file [system format] [-option parameter] [disk name]`
File system formats are: mkfs.cramfs, mkfs.ext2, mkfs.ext3, mkfs.msdos, mkfs.vfat

6. Disk verification: fsck command
Command like: `fsck -t file [system format] [-option parameter] [disk name]`
File system formats are: fsck.cramfs, fsck.ext2, fsck.ext3, fsck.msdos, fsck.vfat

## Bioinformatics on linux cases

#### 1: install software

1. Download software:
`wget –c ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz`

2. Unzip file: `tar zxvf ncbi-blast-2.6.0+-x64-linux.tar.gz`

3. Install the software:
  `cd soft`
  `./configure`
  `make`
  `make install`
    
4. Update PATH: 
  Add the `./ncbi-blast-2.6.0+/bin` directory to the environment variables: `vim ~/.bashrc`
  Add the following statement, save and exit: `export PATH=./ncbi-blast-2.6.0+/bin:$PATH`
  Update environment variables: `source ~/.bashrc`

#### 2: download data from database

1. wget command
  Command like: `wget [-option parameter] [URL]`
  - `-c` When connection has been cut off,you can use this parameter to resume the transfer.
  - `-i` When there are multiple files to download, you can write the URL of each file in a download.txt. Then type command:`wget -i download.txt`.
  - `-r` Download all files from this website including all addresses pointed to by the site.
  For example: `wget [https://www.ncbi.nlm.nih.gov/sra/SRR8467693](https://www.ncbi.nlm.nih.gov/sra/SRR8467693)`
2. sratoolkit
  1. Install sratoolkit
  2. Write command like:
  `module load sratoolkit
  prefetch SRR8467686`

#### 3: fastq-dump

Command like: `fastq-dump [-option parameter]`
- `--split-spot`: Split paired-end sequencing into two copies, but put them in the same file
- `--split-files`: Divide paired-end sequencing into two copies and put them in different files, but discard the reads that one party has but one does not.
- `--split-3`: Divide the paired-end sequencing into two copies and put them in different files, but the reads that one party has and the other does not will be put in a separate folder
- `-o` Output path