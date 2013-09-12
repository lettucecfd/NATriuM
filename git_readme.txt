This readme describes the creation and working with a GIT-Repository

////////////////////////////////
// CREATION OF GIT REPOSITORY //
////////////////////////////////

SET UP SERVER
1) Enable ssh login without password
 - stay on local machine
 - check if local host has rsa key (if yes, there are files ~/.ssh/id_rsa and ~/.ssh/id_rsa.pub)
   ~/.ssh/id_rsa contains the SECRET private key -> NEVER USE SOMEWHERE ELSE THEN ON LOCAL MACHINE;
   ~/.ssh/id_rsa.pub is the public key
 - if not existing, create rsa-key on local host (ssh-keygen -t rsa)
 - Login on server (e.g. ssh akraem3m@wr0.wr.inf.h-brs.de)
 - create a file ~/.ssh/authorized_keys2 (if not exists)
 - add the public key of the local host to ~/.ssh/authorized_keys2 on the server
   (if that does not work later, try to append to ~/.ssh/authorized_keys)
 - make sure that the following permissions are set (on the server):
   1) only you have write access to your home
   2) chmod 640 ~/.ssh/authorized_keys
   3) chmod 700 ~/.ssh
 - If everything went right, you don't have to enter a password anymore when logging in to the server
 
2) Create git hub
 - create directory (e.g. mkdir /scratch/akraemer/NATriuM.git)
 - enter directory 
 - Initialize git (git init --bare)
   If .git file already existed in the directory, you have to delete it (rm -rf .git)
 - git update-server-info # If planning to serve via HTTP
 
3) Initial commit of local files
 - Install git on local host (sudo apt-get install git)
 - make sure the name of the directory is the same as on the repository (e.g. NATriuM)
 - enter directory
 - Initialize git on local host (git init)
 - Add all files of the directory to repository (git add *)
 - Commit (git commit -m "Initial commit message")
 - Add remote repository as origin (git remote add origin user@server:/path/to/repository)
   e.g. git remote add origin akraem3m@wr0.wr.inf.h-brs.de:/scratch/akraemer/NATriuM.git
 - Push the files to repository (git push -u origin master)

4) Make project "shared" in eclipse
 - Right click on project -> Team -> Make shared -> ...

/////////////////////////////////
// WORKING WITH GIT REPOSITORY //
/////////////////////////////////
1) Upgrade to newest eclipse release
 - Start eclipse as root (e.g. gksu /opt/eclipse/eclipse)
 - Make shure the Proxy settings are right (Window->Preferences->General->Network Connections)
   (No proxy server for socks protocol!!!)
 - Add current release's site to updates (Window->Preferences->Install/Update->Add)
   (e.g. http://download.eclipse.org/releases/kepler)
 - Help->Check for updates... -> Install the new version

2) Add EGit to ecplise
 - Add egit site to updates (Window->Preferences->Install/Update->Add)
   (http://download.eclipse.org/egit/update)
 - Help->Install new software
   Work with ->  http://download.eclipse.org/egit/update
   Install Eclipse Git Team Provider
   
3) Configure EGit
 - Window->Preferences->Team->Git->Configuration
 - Add.. key: user.name value: Your Name
 - Add.. key: user.email value: Your Email
 - Define git repository (Git->Default repository folder)
   (e.g. /home/kraemer/git)
 - Create .gitignore-file in Project folder to specify files and directories that will not be synced
   Use wildcards; or whole directories 

 
 

/////////////////////////////////
// WORKING WITH GIT REPOSITORY //
/////////////////////////////////

ADD USER


