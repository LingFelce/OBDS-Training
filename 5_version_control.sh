# Starting out
# navigate to appropriate folder

mkdir - p devel/example-repo # make new repository
cd devel/example-repo # move to repository
ls -al # check files
git status #check repository status

git init # initialise repository
ls - al # should now see .git folder
git status # check git status
git log # check git log
git diff # check file history
git add # staging
git commit # commit with comment

# git remote - put your local repository on github
# create your new repository on github.com - DO NOT INITIALISE WITH README FILE!
git remote -v # see remotes
git remote add origin git@github.com:<user-name>/<name_of_repository>.git #add remote
git push -u origin master # put code on github
# make changes online
git pull origin master # get code from github

# git clone - setting up a copy of an online repository to your computer
git clone git@github.com:<link_here>
git remote -v
git push origin master
git pull origin master

# Notes on using communal repository:
git pull # to get latest version
git branch <name> # create new branch of local repository - do all work on this branch
git checkout <name> # switch to new branch
git add
git commit -m "update comment"
git push origin <name> # push branch to git hub
# on github website, select branch and select new pull request. Merge branches and delete branch.
git checkout master # switch to master branch
git branch -d <name> # delete new branch
git pull # pull latest version to local master repository
