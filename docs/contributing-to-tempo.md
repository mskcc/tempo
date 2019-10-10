# Contributing to Tempo

This document provides a brief set of guidelines for contributing to Tempo.\
<small>We borrowed heavily from the guidelines used by the [cBioPortal](https://github.com/cBioPortal/cbioportal) for third-party contributions, as outlined [here](https://github.com/cBioPortal/cbioportal/blob/master/CONTRIBUTING.md).</small>

### Background
We use the "fork and pull" model for collaborative software development, as explained in the [GitHub Help Page of Using Pull Requests](https://help.github.com/articles/using-pull-requests/):

>_The fork & pull model lets anyone fork an existing repository and push changes to their personal fork without requiring access be granted to the source repository. The changes must then be pulled into the source repository by the project maintainer. This model reduces the amount of friction for new contributors and is popular with open source projects because it allows people to work independently without upfront coordination._

Users are welcomed and encouraged to submit changes they would like to see to Tempo!

### Contributing Code Changes: Making a Pull Request

Once you have forked the repo, you need to create your code contributions within a new branch of your forked repo. For general background on creating and managing branches within GitHub, see: [Git Branching and Merging](https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging).

* To begin, create a topic branch from where you want to base your work.
* For a new feature, this is usually the branch **develop**.  We regularly release changes in **develop** into a formal versioned release in the **master** branch.

Users usually create a branch like so:

```
git checkout develop
git checkout -b [name_of_your_new_branch]
```

Then commit any code changes and push your branch back to GitHub like this:

```
git add . 
git commit -m "give a detailed commit message here"
git push origin [name_of_your_new_branch]
```

A few tips:

* Please name your branch and all commits with descriptive names describing your changes.
* Make commits in logical/cohesive units.
* Make sure you have added the necessary tests for your changes.
* Please make small pull requests for review. 

When you are ready to submit your pull request:

* Push your branch to your GitHub project.
* Open a pull request on GitHub to the **develop** branch at mskcc/tempo for a new feature.

We will then review your suggested changes, and possibly integrate these. 

For instructions on submitting a pull request, please see: [Using Pull Requests ](https://help.github.com/articles/using-pull-requests/) and [Sending Pull Requests](http://help.github.com/send-pull-requests/).

