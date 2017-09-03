# Contribute to ChemApps

ChemApps provides a suite of utility tools and programs for thermochemistry
and chemical kinetics. This document describes how you can contribute.

## How can I contribute?

### Reporting bugs

You can report bugs by [open up a new issue](https://gitlab.com/stigrs/chemapps/issues/new) 
on GitLab. Please ensure that the bug not already reported by searching under
[Issues](https://gitlab.com/stigrs/chemapps/issues). Be sure to include a *title 
and clear description*, as much relevant information as possible, and a *code 
sample* or an *executable test case* demonstrating the expected behavior that is 
not occurring. 

### Fixing bugs

Open a new GitLab pull request with the patch. Ensure that the pull request
description clearly describes the problem and the solution. Include the relevant
issue number if applicable. Before submitting, please read the coding style 
guide to know more about coding conventions.

### Suggesting enhancements

Suggest your change and start writing code in accordance with the coding 
style guide if you get positive feedback from the main developer @stigrs. 
**Never commit any changes to the master or the release branches!** Instead,
do this:
1. Create a fork
2. Clone the fork
3. Create a feature branch
4. Commit to the feature branch
5. Push the feature branch to the fork
6. Submit a merge request, making sure that:
    * Unfinished work is not submitted
    * The commits are atomic and the commit messages are useful
    * The new code is tested and validated
    * The new code is documented

## Coding Style Guide
* The source code should follow *C++ Core Guidelines*;
see http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines for details. 
* The source code shall follow the style supplied in the .clang-format file.
* New features shall be tested and validated with [Catch](https://https://github.com/philsquared/catch)
test cases.
