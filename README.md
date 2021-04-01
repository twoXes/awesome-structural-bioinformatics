# Awesome Protein Folding

![chemistry](chemistry.gif)

Want a better way to learn than some random repo on github? 

Spend 4-12 years of your life and hundreds of thousands of dollars chasing a paper with a stamp on it 🥇.

Or feed yourself 🍼.

## Why Awesome Protein Folding? 

Books on Cheminformatics, Bioinformatics, Quantum Chemistry strangle the subject to sleep 😴 and command a wild price 🤑 for the naps they induce. 

Knoweldge should be cheap, fast enjoyable, silly, shared, disproven, contested, and most of all free. This book combats knowledge hodlers, and innovation stifflers by not being boring and old. This is for the young of mind and young of spirit 🚼 that love to fold. 

## Other Free Books You Should Read Instead of This Repo

[Chemisty 2E](https://openstax.org/details/books/chemistry-2e) - :atom: Equivalent to 201 & 202 Level Chemistry Book

[Chemistry: Atoms First 2E](https://openstax.org/details/books/chemistry-atoms-first-2e) :atom: Fork of 2E but not with more Atoms!!!!

[Biology 2E](https://openstax.org/details/books/biology-2e) 👽 Like Chemistry 2E but Biology 

[Artificial Intelligence: A Modern Approach]( https://github.com/aimacode/aima-python) 🤖 The Gospel of Machine Learning

[Neural Networks and Deep Learning](http://neuralnetworksanddeeplearning.com/) 🤖 Michael Nielsen writes another masterpiece - About Deep Learning - if you are into that sort of thing. 

[Reinforcement Learning](http://incompleteideas.net/book/RLbook2020.pdf) 🤖  The only book you need on the subject 

[Pattern Recognition and Machine Learning](https://www.microsoft.com/en-us/research/uploads/prod/2006/01/Bishop-Pattern-Recognition-and-Machine-Learning-2006.pdf) 🤖 Another classic banger

## Why Should You Care

![Prometheus](prometheus.gif)

At their core, Crick and thus *Prometheus* are taking about evolution at the molecular level. 

 > Structure-function relationships are the fundamental object of knowledge in protein chemistry; they allow us to rationally design drugs, engineer proteins with new functions, and understand why mutations cause disease. [source](https://www.chemistryworld.com/features/on-the-origin-of-proteins/3004719.article)

## Scoring Function 

![Goal](goal.gif)

 There are two common approaches to building a score function: 
* **potentials of mean force**
    * often called statistics- or Boltzmann-based force fields
    * measuring distance as a reflection of statistical tendencies within proteins
    * . One takes a large set of proteins, collects statistics and converts them to a score function. One then expects this function to work well for proteins not included in its parameterisation. 
* **an optimization calculation**
   * select underlying basis function 
      * quasi-Lennard-Jones 
      * various sigmoidal functions 
   * We can say that the correct structure is whatever is given in the protein data bank, but unfortunately, there is almost an infinity of incorrect structures for a sequence and one would like the score function to penalize all of them
  *  One way to encode this idea is to adopt a statistical approach and try to consider the distribution of incorrect structures
 [source](http://web.stanford.edu/class/cs273/refs/torda_chapter_proteomics.pdf)

Allowing gaps and insertions at any position and of any length leads to a combinatorial explosion of possibilities. The calculation can be made tractable by restricting the search space and forbidding gaps except in recognised loops in template structures.

There is a score function and a fast method for producing the best possible sequence to structure alignments and thus the best models possible. Unfortunately, the problem is still not solved

### Probability 

![Probably](probably.gif)

### Tools 

[OpenMM](http://openmm.org/) - molecular simulation


(2020) [Deep Learning Methods in Protein Structure Prediction](https://www.sciencedirect.com/science/article/pii/S2001037019304441)
