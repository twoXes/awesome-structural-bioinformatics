import streamlit as st
import Bio
from Bio.Blast import NCBIWWW

'''

#  Molecular Informatics Artifically Deep

An artifically deep dive into machine learning methods for molecular informatics 

'''

st.image("chemistry.gif")

'''
Most code written in Python :dragon:  except when it isn't. 

Want a better way to learn? Spend 4-12 years of your life and hundreds of thousands of dollars chasing a paper with a stamp on it :medal:

Or feed yourself :baby_bottle:

## Why This Book

Most books on Cheminformatics, Bioinformatics, Quantum Chemistry, are good for inducing sleep 😴 and  fleecing your wallet :moneybag: 

Knoweldge should be cheap, fast enjoyable, silly, shared, disproven, contested, fun and most of all free. This book combats knowledge hodlers, and innovation stifflers by not being boring and old. 

Though we look at  supposedly outdated (if you can't date yourself, who can you date?) Machine Learning methods as applied to Molecular Informatics, this is for the young of mind and young of spirit 🚼 

## Other Free Books You Should Read Instead of This

[Chemisty 2E](https://openstax.org/details/books/chemistry-2e) - :atom_symbol: Equivalent to 201 & 202 Level Chemistry Book

[Chemistry: Atoms First 2E](https://openstax.org/details/books/chemistry-atoms-first-2e) :atom_symbol: Fork of 2E but not with more Atoms!!!!

[Biology 2E](https://openstax.org/details/books/biology-2e) 👽 Like Chemistry 2E but Biology 

[Artificial Intelligence: A Modern Approach]( https://github.com/aimacode/aima-python) 🤖 The Gospel of Machine Learning

[Neural Networks and Deep Learning](http://neuralnetworksanddeeplearning.com/) 🤖 Michael Nielsen writes another masterpiece - About Deep Learning - if you are into that sort of thing. 

[Reinforcement Learning](http://incompleteideas.net/book/RLbook2020.pdf) 🤖  The only book you need on the subject 

[Pattern Recognition and Machine Learning](https://www.microsoft.com/en-us/research/uploads/prod/2006/01/Bishop-Pattern-Recognition-and-Machine-Learning-2006.pdf) 🤖 Another classic banger

## YouTube University

Don't wait for "school," Attend YouTube U for any prereqs you need to get started, today:

What is the the structure of DNA?

'''

st.video('https://youtu.be/o_-6JXLYS-k')
st.video('https://youtu.be/MODnIkQvyz0')
'''
## The Rise of Data

1984 -  [Protein Information Resource](https://proteininformationresource.org/)

12

Structural Classification of Proteins [Database](http://scop.mrc-lmb.cam.ac.uk/)

## [Protein Folding](https://sitn.hms.harvard.edu/flash/2010/issue65/)
'''

st.video('https://youtu.be/hok2hyED9go')

st.video('https://youtu.be/1peFJ_-N7V8')

'''
The [three questions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2443096/) of the Protein Folding problem:
* What is the folding code?
* What is the folding mechanism?
* Can we predict the native structure of a protein from its amino acid sequence?

Protein Structure:
 * sequence of amino acids
 * secondary structure
    * a-Helixs
    * B-sheets
 * tertiary structure
 * Quaternary structure
    * subunit organization

The biggest factor in a proteins ability to fold is the thermodynamics of the structure.[source](https://chem.libretexts.org/Bookshelves/Biological_Chemistry/Supplemental_Modules_(Biological_Chemistry)/Proteins/Protein_Structure/Protein_Folding)
 * if a proteing is not in the lowest energy conformation it will continue to move and adjust until it finds its most stable state.

Another factor is hydrophic interactions within the protein
 * hydrophobic interactions have an impact not just on the primary structure but then lead to changes seen in the secondary and tertiary structure as well.
 * The hydrophobic interactions are shown to have an impact on the protein even after it has found the most stable conformation in how the proteins can interact with each other as well as folding themselves.

Another type of interaction seen when the protein is folding is the disulfide linkages that form in the protein a  sulfur- sulfur chemical bond that results from an oxidative process that links nonadjacent (in most cases) cysteine’s of a protein.

### AlphaFold
'''
st.image('alphafold.gif')


