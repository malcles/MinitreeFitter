MinitreeFitter
==============


Saclay group MinitreeFitter: make workspace and dataset for Hgg analysis

Install
--------

** replace username in the following by your actual git user name **

1. Fork MinitreeFitter in https://github.com/fcouderc/MinitreeFitter
2. Go to an lxplus machine (SLC6 but SLC5 should work just as good)
3. Download MinitreeFitter:
- git clone git@github.com:fcouderc/MiniTreeFitter MiniTreeFitter
- cd MiniTreeFitter
- source etc/scripts/setup.sh
- make -j 4   
4. add your remote fork working directory
- git remote add minitreeFork https://github.com/cms-analysis/flashgg
5. some help
- MiniTreeFitter --help


Committing 
----------

git is managing everything locally, once you are happy with a snapshot of your work, push it to github so others can use it

1. adding and committing
- git add <files> <dirs>   /// files or dirs that have been changed and that you want to commit
- git commit -m "message"  /// commit what was added
- git status               /// check the status of the snapshot
- git remote -v            /// gets the info of the remote github dir, useful to get the ALIAS of these directories
       	                    /// (the first name is probably origin, second one is probably minitreeFork)
   git push ALIAS BRANCH (to push your snapshot, BRANCH probably master)
2. fetching from remote dir
- git fetch ALIAS BRANCH   
3. ask for a pull request:
- push snapshot in the fork aliased minitreeFork
  git push minitreeFork master
- ask for pull request on the web

Working example
---------------

** always do setup before anything **
- source etc/scripts/setup.sh
- MiniTreeFitter -d /afs/cern.ch/work/f/fcouderc/public/minitreeFitterExample/cic


