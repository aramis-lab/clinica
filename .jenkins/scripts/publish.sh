#!/bin/sh
# Script written by Thierry Martinez and modified to be used in this
# repository.
set -ex
# The script iterates each time on all branches to compile the pages.
# The artifact needs to contain all the pages and I don't see any
# mean to combine all the artifacts produced by all the branches 
# (There is the "needs:" field in .gitlab.yml to invoke an artifact
# produced by an other branch, but this would require to hard-wire all
# the branch names in .gitlab.yml...)

branch=$1
mkdocs build || true
mv site "$branch"

## Remotes refs are of the form "refs/remotes/origin/branch-name"
## so we strip the three first path components to keep only "branch-name".
## We don't keep neither "dont_publish.*" branches nor HEAD branch.
#branches="`\
#  git branch --format='%(refname:lstrip=2)' |\
#  grep -v '^dont_publish\|^HEAD$'`"
## The branches/ directory will temporarily contain one directory with
## the generated pages for each branch
#mkdir branches
#for branch in $branches; do
#  # We are in detached HEAD: use "git reset" to force getting the branch
#  git reset --hard origin/$branch
#  mkdocs build || true
#  mv site branches/$branch
#done
## The master branch becomes the public root, and other branches go in
## the branches/ subdirectory of public.
#mv branches/dev public
#mv branches public/

