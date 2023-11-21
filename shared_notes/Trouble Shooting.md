>Created by Zitian.
>Last modified [11-21-2023](11-21-2023).


#### **Purpose:** this document is a collection of all errors that were encountered while running [WES pipeline](./RoTATion%20-%20Germline%20Varient%20Calling%20WES.md) AND has been reported to Zitian (:D). Therefore, this doc contains solutions to some commonly-seen errors (but not all!).
---


## Error Type 1: File Not Found

The first type of error you'll see very often is that certain files could not be found by your script. Error message is usually like this:
```bash
"A USER ERROR has occurred: Couldn't read file [filename]. Error was: It doesn't exist."

"Fail to open file [filename]: No such file or directory."

"FileNotFoundError: [Errno 2] No such file or directory: [filename]."
...
```

## Solution to Error Type 1

Check the following stuff **before** you reach out for help:

### a. Are you in the correct directory?

`ls` in side your terminal/console window, see whether the file `[filename]` do exist. If not, make sure to `cd` to the correct directory before executing any jobs.

### b. Have you exported LSF docker volumes?

Export the LSF docker volumes **before** you request a docker ensures whatever you're going to execute (code, scripts, etc.) _sees the files inside the directories that have been exported_. 

For example, if we're only exporting `jin810/Active/testing`, then the code will not be able to see any files inside `jin810/Active/Projects`. In this case, if you're, say, using anything in `Projects` folder as inputs to any of your scripts, you will very likely get a `file not found` error, as RIS is defaulted to be an empty environment and needs you to tell it what to do and in what space to search files.

**_You can export LSF docker volumes by running:_**
```bash
export LSF_DOCKER_VOLUMES="/storage1/fs1/bga/Active:/storage1/fs1/bga/Active /scratch1/fs1/jin810:/scratch1/fs1/jin810 /storage1/fs1/jin810/Active:/storage1/fs1/jin810/Active $HOME:$HOME"
```

### c. Is your file executable?

By default, all scripts are _non-executable_ upon creation. If you're trying to execute such scripts, the RIS will throw an `file not found` error on you. 

**_To make them executable, you can run:_**
```bash
chmod +x [filename].[file_extension]
# for example: chmod +x vcf2bed.py will make this python file executable.
```

### d. Do you have the correct permission?

Each file has designated permissions (read, write, and execute) for the users (you), the group (others in compute-jin810), and others (other RIS users). 

Most commonly when you **copied a script from other people's directory**, you will likely NOT have a permission to write or execute the file, which will give you a similar error as described in [the above section](#c-is-your-file-executable).

Two potential solutions here. You can either create a new file and copy the code over; or change the permissions to the file using various `chmod` commands. One commonly used one is `chmod 755`, which grants you read, write, and execute permissions to the file.


## Error Type 2: Package Not Found

As different people configured their local and/or RIS environment differently, there is a chance that scripts successfully ran on other people's machine throw errors at you. It is totally common, so don't be over frustrated!

The error message will most likely be something like:
```bash
"bash: [package_name]: command not found"

"/some/path: cannot find [package_name]"
...
```

## Solution to Error Type 2

This issue is most likely associated with your environmental variable, `PATH`. You can see what directories is in your `PATH` by doing `echo $PATH` in your terminal. This variable defines under which directories the docker you're going to request will search for packages.

This is usually not a problem if you've built your own docker, but in many cases we are using dockers built by someone else. 

```bash
## Request the docker you will be using in interactive mode:
bsub -Is -G compute-jin810 -q general-interactive -a 'docker([docker_name])' /bin/bash

## Search for paths to packages
whereis [package_name]    # find places contain the package and its dependencies
which [package_name]      # find the executable path to the package

## Exit interactive mode or open a new terminal window/tab
export PATH="/path/of/package:$PATH"

## Then submit jobs as you would usually do
bsub ...
```

For example, if I got an error message saying `python not found`, and through `which python`, I figured that python is under `opt/conda/bin/`. Then, I can do `export PATH="/opt/miniconda/bin:$PATH"` before I run `bsub` in order for the script to be able to find python.


## Error Type 3: Disk Quota Exceeds

Error message like this:
```bash
"write /dev/stdout: disk quota exceeded"
```

## Solution to Error Type 3

Check your home directory and see if it is close to **10GB**.

```bash
du -sh ~            # this will give you the entire size of your home directory

# if the above command gives a result close to 10GB:

du -sh $(ls -A)     # this will give you the size of each directory under your $HOME, including hidden ones
```

Then, find the directories with largest sizes (usually hidden directories), and ask anyone (preferably someone working on computational stuff) and/or google on whether it is safe to delete that folder. 

