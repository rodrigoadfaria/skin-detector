# Combined Rules for Skin Detection
This repo is host for the experiments with combined rules to detect skin color pixels in colored images.

## How do I get setup?
You can open this project within Visual Studio or compile the source code using gcc in Linux/Mac

## What is in each branch?
We have four different branches, each one with:

* master
  
  The original source code from Nadia Brancati paper.

* reverse
  
  The source code changed to refactor the correlation rules to put them in terms of the estimated value of Pcr.

* combined

  The combination of the rules given in the original (master branch) and reverse hipothesis.

* combined_gs

  This branch we are combining different pairs of Y min/max percentils in order to find out the best combination of them.

* neighbors
  
  A neighborhood extended method based on the combined rules.

## Who do I talk to?
### Repo owner or admin

Rodrigo Faria, rodrigoadfaria@gmail.com
