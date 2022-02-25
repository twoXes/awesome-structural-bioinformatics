# Awesome Structural BioInformatics [![Awesome](https://cdn.rawgit.com/sindresorhus/awesome/d7305f38d29fed78fa85652e3a63e154dd8e8829/media/badge.svg)](https://github.com/sindresorhus/awesome)

*A curated list of awesome structural bioinformatics frameworks, libraries, software and resources.*

> So let it not look strange if I claim that it is much easier to explain the movement of the giant celestial bodies than to interpret in mechanical terms the origination of just a single caterpillar or a tiny grass. - Immanuel Kant, *Natural History and the Theory of Heaven*, 1755

Books on Cheminformatics, [Bioinformatics](http://thegrantlab.org/teaching/material/Intro_Structural_Bioinformatics_Barry.pdf), Quantum Chemistry strangle the subject to sleep 😴 and command a wild price 🤑 for the naps they induce.

Want a better way to learn than some random repo on github? 

Spend 4-12 years of your life and hundreds of thousands of dollars chasing a paper with a stamp on it 🥇.

Or feed yourself 🍼.

Information should be cheap, fast enjoyable, silly, shared, disproven, contested, and most of all free.

Knowledge hodlers, and innovation stifflers are boring and old. This is for the young of mind and young of spirit 🚼 that love to dock & fold. 

[Proteins[(https://github.com/Ramb0a/awesome-structural-bioinformatics#protein-bioinformatics)

[Genetics](

# Protein BioInformatics

## Protein Folding

 > Structure-function relationships are the fundamental object of knowledge in protein chemistry; they allow us to rationally design drugs, engineer proteins with new functions, and understand why mutations cause disease. [- On The Origin of Proteins](https://www.chemistryworld.com/features/on-the-origin-of-proteins/3004719.article)

> There is now a testable explanation for how a protein can fold so quickly: A protein solves its large global optimization problem as a series of smaller local optimization problems, growing and assembling the native structure from peptide fragments, local structures first. [- The Protein Folding Problem](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2443096/)

![Simpsons](simpsons.gif)

The protein folding problem consists of three closely related puzzles:
* (a) What is the folding code? 
* (b) What is the folding mechanism?
* (c) Can we predict the native structure of a protein from its amino acid sequence? [source](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2443096/)


### Data Sources

![Prometheus](prometheus.gif)

[CATH/Gene3D](https://www.cathdb.info/) - 151 Million Protein Domains Classified into 5,481 Superfamilies

[NCBI Conserved Domains Database](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi) - resource for the annotation of functional units in proteins

[Protein Data Bank](https://www.rcsb.org/)

[Scop 2](https://scop.mrc-lmb.cam.ac.uk/) - Structural Classification of Proteins

[UniProt](https://www.uniprot.org/) -  comprehensive, high-quality and freely accessible resource of protein sequence and functional information.

[Fold@Home](https://foldingathome.org/about/)

### Deep Learning Protein Folding

#### [AlphaFold 14](https://www.predictioncenter.org/casp14/doc/presentations/2020_12_01_TS_predictor_AlphaFold2.pdf)

  * [:computer: Code](https://github.com/deepmind/alphafold)
  * [:book: Paper 2](https://www.nature.com/articles/s41586-021-03819-2_reference.pdf)
  * [:book: Paper](https://www.nature.com/articles/s41586-019-1923-7.epdf?author_access_token=Z_KaZKDqtKzbE7Wd5HtwI9RgN0jAjWel9jnR3ZoTv0MCcgAwHMgRx9mvLjNQdB2TlQQaa7l420UCtGo8vYQ39gg8lFWR9mAZtvsN_1PrccXfIbc6e-tGSgazNL_XdtQzn1PHfy21qdcxV7Pw-k3htw%3D%3D)
  * [:newspaper: article](https://deepmind.com/blog/article/AlphaFold-Using-AI-for-scientific-discovery)
  * [AlpahFold 14 Results Discussion](https://dasher.wustl.edu/bio5357/discussion/oxford-alphafold2.pdf)
  * [What AlphaFold means for Structural BioInformatics](https://ammiellewb.medium.com/what-alphafold-means-for-structural-bioinformatics-78117adb7d11)
  * [AlphaFold 2 Explained](https://youtu.be/B9PL__gVxLI) - Yanick Video
  * [Illustrated Transformer](kjalammar.github.io/illustrated-transformer/)
  * [Transformers from Scratch](http://peterbloem.nl/blog/transformers)

#### [AlphaFold 13](https://www.predictioncenter.org/CASP13/doc/presentations/Pred_CASP13-Structure-AlphaFold-Jumper.pdf)

  * [:floppy_disk: Code](https://github.com/deepmind/deepmind-research/tree/master/alphafold_casp13)
  * [:floppy_disk: Code](https://github.com/dellacortelab/prospr) - Prospr - Open Source Implementation
  * [:book: Prospr Paper](https://www.biorxiv.org/content/10.1101/830273v1) 
  * [AlphaFold @ Casp13: What Just Happened?](https://moalquraishi.wordpress.com/2018/12/09/alphafold-casp13-what-just-happened/) 

#### [MiniFold](https://github.com/hypnopump/MiniFold) - Open Source toy example of AlphaFold 13 algorithm 

> The DeepMind work presented @ CASP was not a technological breakthrough (they did not invent any new type of AI) but an engineering one: they applied well-known AI algorithms to a problem along with lots of data and computing power and found a great solution through model design, feature engineering, model ensembling and so on...

> Based on the premise exposed before, the aim of this project is to build a model suitable for protein 3D structure prediction inspired by AlphaFold and many other AI solutions that may appear and achieve SOTA results.

![MiniFold](minifold.png)

> Two different residual neural networks (ResNets) are used to predict angles between adjacent aminoacids (AAs) and distance between every pair of AAs of a protein. For distance prediction a 2D Resnet was used while for angles prediction a 1D Resnet was used.

#### PDNet

> As deep learning algorithms drive the progress in protein structure prediction, a lot remains to be studied at this merging superhighway of deep learning and protein structure prediction. Recent findings show that inter-residue distance prediction, a more granular version of the well-known contact prediction problem, is a key to predicting accurate models. However, deep learning methods that predict these distances are still in the early stages of their development. To advance these methods and develop other novel methods, a need exists for a small and representative dataset packaged for faster development and testing. In this work, we introduce protein distance net (PDNET), a framework that consists of one such representative dataset along with the scripts for training and testing deep learning methods. The framework also includes all the scripts that were used to curate the dataset, and generate the input features and distance maps.

[ 💻 Github](https://github.com/ba-lab/pdnet/)

[ :book: Paper](https://www.nature.com/articles/s41598-020-70181-0) 

[ :vhs: YouTube](https://youtu.be/uAIuA1O7iE8)


## Protein - Ligand Docking


"Docking is a method which predicts the prefered orientation of one molecule to a second when bound to each other to form a stable complex. Knoweldge of the prefered orientation in turn may be used to predict the strength of association or binding affinity between two molecules using scoring functions."

* Pose - A conformation of the receptor and ligand molecules showing some intermolecular interactions (which may include hydrogen bonds as well as hydrophobic contacts
* Posings - The process of searching for a pose in which there are favorable interactions between the receptor and the ligand molecules.
* Scoring - The process of evaluating a particular pose using a number of descriptive features like number of intermolecular interactions including hydrogen bonds and hydrophobic contacts.

* The best docking algorithm should be the one with the best scoring function and
the best searching algorithm
[source](https://www.beilstein-journals.org/bjoc/articles/12/267)
* No single docking methods performs well for all targets and the quality of
    docking results is highly dependent on the ligand and binding site of
    interest [source](https://www.beilstein-journals.org/bjoc/articles/12/267)


In the early 1990s many approved HIV protease inhibitors were developed to target HIV infections using structure-based molecular docking. [source](https://www.beilstein-journals.org/bjoc/articles/12/267)
* Saquinavir
* Amprenavir 

Scoring Functions in MD can be categorized into:
* **knowledge based** - stastical potentials, frequency of interaction occurance,
    Boltzmann distribution, dataset dependent  
* **force-field based** - energy functions via molecular mechanics, coulombic
    interactions, van der Waals interactions (Lennard-Jones potential)
      * CHARMM (chemistry at Harvard macromolecular mechanics)
      * AMBER (assisted model building and energy refinement)
* **empirical** - binding free energy calculated as the weighted sum of
    unccorrelated terms,(example - hydrogen bonds, hydrophobicity), Regression
    analysis find the best weights for each term
      * HYDE (part of BioSolveIT tools)
      * ChemScore
      * SCORE
* **consensus** - combines scoring functions types into ensemble 
     
      * X-CSCORE
      * MultiScore



One of the first appearances of Molecular Docking is said to have been 1982's 

### [A Geometric Approach to MacroMolecule Ligand Interactions](https://www.sciencedirect.com/science/article/abs/pii/002228368290153X#:~:text=A%20geometric%20approach%20to%20macromolecule%2Dligand%20interactions%E2%98%86&text=Algorithms%20are%20presented%20that%20examine,sites%20on%20a%20macromolecular%20surface.)

They tell us Molecular Docking = "To position two molecules so that they interact favorably with one another..."

How???

> Our approach is to reduce the number of degrees of freedom using simplifying
assumptions that still retain some correspondence to a situation of biochemical
interest. Specifically, we treat the geometric (hard sphere) interactions of two rigid
bodies, where one body (the “receptor”) contains “pockets” or “grooves” that form
binding sites for the second object, which we will call the “ligand”. Our goal is to fix
the six degrees of freedom (3 translations and 3 orientations) that determine the
best relative positions of the two objects.

> Does the program reproduce known ligand-receptor geometries? If so, does it also provide alternative structures that are geometrically
reasonable? To these ends, we have examined two systems for which the ligand receptor
geometry has been established by crystallographic means.

What is the result of this Docking? 

> (1) Structures quite near the “correct” structures are readily recovered and
identified as feasible solutions.
> (2) Other families of structures are found that are geometrically reasonable and
that can be tested by simple scoring schemes, chemical intuition, or visual
inspection with computer graphics.

> Without allowing molecular flexibility, many aspects of ligand-receptor interactions are not properly described.

---

> A common approach to docking combines a scoring function with an optimization algorithm. The scoring function quantifies the favorability of the protein-ligand interactions in a single pose, whichcan be conceptualized as a point in a continuous conformation space. A stochastic global optimization algorithm is used to explore and sample this conformation space. Then, local optimization is employed on the sampled points, usually by iteratively adjusting the pose in search of a local extremum of the scoring function. Ideally, the scoring function is differentiable to support efficient gradient-based optimization.

> The information obtained from the docking technique can be used to suggest the binding energy, free energy and stability of complexes. At present, docking technique is utilized to predict the tentative binding parameters of ligand-receptor complex beforehand.

> There are various databases available, which offer information on small ligand molecules such as CSD (Cambridge Structural Database), ACD (Available Chemical Directory), MDDR (MDL Drug Data Report) and NCI (National Cancer Institute Database).


### Scoring Function 

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

## Protein - Ligand Docking Tools

*Tools for exploring how two or more molecular structures fit together*

[AutoDock](http://autodock.scripps.edu/) - suite of automated docking tools designed to predict how small molecules bind to a receptor of known 3D structure

[AutoDock Vina](http://vina.scripps.edu/) - significantly improves the average accuracy of the binding mode predictions compared to AutoDock

* [:book: Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/)

[Gnina](https://github.com/gnina/gnina) - deep learning framework for molecular docking -inside deepchem (/dock/pose_generation.py)

[GOMoDo](https://gomodo.grs.kfa-juelich.de/php/about.php) - GPCR online modeling and docking server

[Smina](https://github.com/mwojcikowski/smina) used for minimization (local_only) as opposed to of docking, makes Vina much easer to use and 10-20x faster. Docking performance is about the same since partial charge calculation and file i/o isn't such a big part of the performance.


## Fusion Proteins 

[ChimPipe](https://github.com/Chimera-tools/ChimPipe#:~:text=ChimPipe%20is%20a%20computational%20method,junctions%20at%20base%2Dpair%20resolution.) - ChimPipe is a computational method for the detection of novel transcription-induced chimeric transcripts and fusion genes from Illumina Paired-End RNA-seq data. It combines junction spanning and paired-end read information to accurately detect chimeric splice junctions at base-pair resolution.

[DeepNF](https://github.com/VGligorijevic/deepNF) - Deep network fusion for protein function prediction | 📖 [paper](https://academic.oup.com/bioinformatics/article/34/22/3873/5026651) 

[DeepPrior](https://github.com/bioinformatics-polito/DEEPrior) - predicts the probability of a gene fusion being a driver of an oncogenic process by directly exploiting the amino acid sequence of the fused protein, and it can prioritize gene fusions from different tumors. Unlike state-of-the-art tools, it also supports easy retraining and re-adaptation of the model | 📖 [paper](https://academic.oup.com/bioinformatics/article/36/10/3248/5722203)

[DeFuse](https://github.com/amcpherson/defuse) - gene fusion discovery using RNA-Seq data. The software uses clusters of discordant paired end alignments to inform a split read alignment analysis for finding fusion boundaries | 📖 [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3098195/)

[FusionCatcher](https://github.com/ndaniel/fusioncatcher) - Finder of somatic fusion-genes in RNA-seq data

[Jaffa](https://github.com/Oshlack/JAFFA) - JAFFA is a multi-step pipeline that takes either raw RNA-Seq reads, or pre-assembled transcripts, then searches for gene fusions

[StarFusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) | 📖 [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1842-9)

# Genomics 

Human genome = 20K+ genes 👖 each responsible for the instructions of building a single protein; encoded within 6 feet of DNA.

DNA 🧬 is the 'instruction manual of life' (thankfully not written by Ikea 🪑 )

Life is encoded in digital form G/C || T/A nucleotide base pairs. 

DNA is encoded ('transcribed) into mRNA and then decoded ('translated) into Proteins. Along the way alot of good and bad stuff happens in the [latent space](https://stats.stackexchange.com/questions/442352/what-is-a-latent-space).

To make the DNA -> RNA transcribing happen - RNA takes a DNA strands - the template strand that runs from 3' to 5' (prime) - and it is uses to form single stranded RNA with all T's replaced by U's. (Thymine replaced by Uracil). The RNA now runs 5' to 3' prime and is identical to the non-template strand DNA sequence (the coding strand), except again Thymine is replaced by Uracil. 

```python

import Bio
from Bio.Seq import Seq

dna = Seq("ACGTTTATCGATCGA")
mRNA = dna.transcribe()
protein = mRNA.translate()

print(dna)
print(mRNA)
print(protein)

>>> ACGTTTATCGATCGA
>>> ACGUUUAUCGAUCGA
>>> TFIDR

```

From the template strand -> RNA strand transcription, translation happens when RNA is turned into proteins via the reading of codons - three letter RNA sequences that encode into specific amino acids. The tranlation looks like this codon table: 


![Codon Table](https://upload.wikimedia.org/wikipedia/commons/thumb/7/70/Aminoacids_table.svg/1920px-Aminoacids_table.svg.png)

Genes come in **Alleles**, the variations of the gene (example gene = hair color, allele = red hair 👩‍🦰 ||  blonde hair 👱‍♂️) 

> Geneomics is the study of all genes in an organism to understand their molecular organization, function, interaction and evolutionary history.

Genomics begins with the discoveries of [Gregor Mendel](https://www.britannica.com/biography/Gregor-Mendel) 

![Griffith Experiment](https://upload.wikimedia.org/wikipedia/commons/thumb/6/6a/Griffith_experiment.svg/1432px-Griffith_experiment.svg.png)


## Genomics Software

[Deep Variant](https://github.com/google/deepvariant) - analysis pipeline that uses a deep neural network to call genetic variants from next-generation DNA sequencing data

[NVIDIA Clara Parabricks Pipelines](https://docs.nvidia.com/clara/parabricks/v3.5/text/software_overview.html) - perform secondary analysis of next generation sequencing (NGS) DNA and RNA data, blazing fast speeds and low cost. Can analyze whole human genomes in about 45 minutes. Includes Deep Variant. 
  * [AWS Clara Parabrick Pipeline Setup](https://aws.amazon.com/marketplace/pp/prodview-apbngojlskcyq) ~ $4 per hour. 
 
 [BioWasm](https://github.com/biowasm/biowasm) - WebAssembly modules for genomics
 
 [FastQ Bio](https://github.com/robertaboukhalil/fastq.bio) - An interactive web tool for quality control of DNA sequencing data
 
 [Minimap2](https://github.com/lh3/minimap2)  - sequence alignment program that aligns DNA or mRNA sequences against a large reference database. For >100bp Illumina short reads, minimap2 is three times as fast as BWA-MEM and Bowtie2, and as accurate on simulated data. | [paper](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)
 
 [gBWT](https://github.com/jltsiren/gbwt) -  graph extension (gPBWT) of the positional Burrows-Wheeler transform (PBWT)
 * [gBWT in Rust](https://github.com/jltsiren/gbwt-rs) 

[VG](https://github.com/vgteam/vg) - tools for working with genome variation graphs

[Cello](https://github.com/CIDARLAB/Cello-v2) - Genetic Circuit Design 

## Genomics Resources

🍼 [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) - develop the technical infrastructure (reference standards, reference methods, and reference data) to enable translation of whole human genome sequencing to clinical practice and innovations in technologies.

[Online Needleman-Wunsch Example](https://berthub.eu/nwunsch/) || [Example II](http://experiments.mostafa.io/public/needleman-wunsch/index.html) || [Great NW Colab](https://colab.research.google.com/github/zaneveld/full_spectrum_bioinformatics/blob/master/content/08_phylogenetic_trees/needleman_wunsch_alignment.ipynb)

# Learning

💭  [Rosalind](https://rosalind.info/problems/locations/)

💭 [Great Introduction to BioInformatics Course](https://gtpb.github.io/ELB19F/) - ELB19F

💭 [Learn BioInformatics in the Browser](https://sandbox.bio/) - Sandbox Bio

☁️ [Biological Modeling](https://biologicalmodeling.org/) - Free Online Course

☁️ [BioInformatic Algorithms Lecture Videos](https://www.bioinformaticsalgorithms.org/lecture-videos)


# Appendix

## Protein References

(2020) [High-Throughput Docking Using Quantum Mechanical Scoring](https://www.frontiersin.org/articles/10.3389/fchem.2020.00246/full)

(2020) [Deep Learning Methods in Protein Structure Prediction](https://www.sciencedirect.com/science/article/pii/S2001037019304441)

(2019)[From Machine Learning to Deep Learning: Advances in scoring functions for protein-ligand docking](https://onlinelibrary.wiley.com/doi/abs/10.1002/wcms.1429)

(2019) [The Unreasonable Effectiveness of Convolutional Neural Networks in Population Genetic Inference, Molecular Biology and Evolution](https://doi.org/10.1093/molbev/msy224)

(2018) [DeepFam: deep learning based alignment-free method for protein family modeling and prediction](https://doi.org/10.1093/bioinformatics/bty275)

(2018) [Derivative-free neural network for optimizing the scoring functions associated with dynamic programming of pairwise-profile alignment](https://doi.org/10.1186/s13015-018-0123-6)

(2017) [Protein-Ligand Scoring with CNN](https://pubs.acs.org/doi/10.1021/acs.jcim.6b00740)

(2017) [Quantum-chemical insights from deep tensor neural networks](https://www.nature.com/articles/ncomms13890)

(2014) [MRFalign: Protein Homology Detection through Alignment of Markov Random Fields](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003500)

(2012) [Molecular Docking: A powerful approach for structure-based drug discovery](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3151162/#:~:text=The%20docking%20process%20involves%20two,assessment%20of%20the%20binding%20affinity.)

(2011) [The structural basis for agonist and partial agonist action on a β(1)-adrenergic receptor](https://pubmed.ncbi.nlm.nih.gov/21228877/)

(2010) [AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/)

(2009) [Amphipol-Assisted in Vitro Folding of G Protein-Coupled Receptors](https://pubs.acs.org/doi/10.1021/bi801729z)

(2005) [GPCR Folding and Maturation](https://link.springer.com/chapter/10.1007/978-1-59259-919-6_3) from [The G Protein-Coupled Receptors Handbook](https://link.springer.com/book/10.1007/978-1-59259-919-6)

(1982) [A Geometric Approach to MacroMolecule Ligand Interactions](https://www.sciencedirect.com/science/article/abs/pii/002228368290153X#:~:text=A%20geometric%20approach%20to%20macromolecule%2Dligand%20interactions%E2%98%86&text=Algorithms%20are%20presented%20that%20examine,sites%20on%20a%20macromolecular%20surface.)

## Genomics References 

(2021) [Re-identification of individuals in genomic datasets using public face images](https://www.science.org/doi/epdf/10.1126/sciadv.abg3296)

(2021) [Accurate, scalable cohort variant calls using DeepVariant and GLnexus](https://academic.oup.com/bioinformatics/article/36/24/5582/6064144)

(2020) [Secure large-scale genome-wide association studies using homomorphic encryption](https://www.pnas.org/content/117/21/11608.long)

(2020) [Optimized homomorphic encryption solution for secure genome-wide association studies](https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-020-0719-9)

(2019) [Are drug targets with genetic support twice as likely to be approved? Revised estimates of the impact of genetic support for drug mechanisms on the probability of drug approval](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008489)

(2018) [A universal SNP and small-indel variant caller using deep neural networks](https://github.com/Ramb0a/awesome-genomics/blob/main/papers/universal_snp_small_indel_variant_caller_using_dnn.pdf)

(2018) [Secure genome-wide association analysis using multiparty computation](https://github.com/Ramb0a/awesome-genomics/blob/main/papers/secure_genome_wide_association_analysis_using_multiparty_computation.pdf)

(2018) [Minimap2: pairwise alignment for nucleotide sequences](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)

(2016) [The sequence of sequencers: the history of sequencing dna](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4727787/)


## Other Free Books You Should Read Instead of This Repo

[Chemisty 2E](https://openstax.org/details/books/chemistry-2e) - :atom: Equivalent to 201 & 202 Level Chemistry Book

[Chemistry: Atoms First 2E](https://openstax.org/details/books/chemistry-atoms-first-2e) :atom: Fork of 2E but not with more Atoms!!!!

[Biology 2E](https://openstax.org/details/books/biology-2e) 👽 Like Chemistry 2E but Biology 

[Artificial Intelligence: A Modern Approach]( https://github.com/aimacode/aima-python) 🤖 The Gospel of Machine Learning

[Neural Networks and Deep Learning](http://neuralnetworksanddeeplearning.com/) 🤖 Michael Nielsen writes another masterpiece - About Deep Learning - if you are into that sort of thing. 

[Reinforcement Learning](http://incompleteideas.net/book/RLbook2020.pdf) 🤖  The only book you need on the subject 

[Pattern Recognition and Machine Learning](https://www.microsoft.com/en-us/research/uploads/prod/2006/01/Bishop-Pattern-Recognition-and-Machine-Learning-2006.pdf) 🤖 Another classic banger

# Brief Explanation of AlphaFold Jax Architecture

AlphaFold2 is Google's state of the art protein structure prediction model.

AF2 predicts 3D coordinates of all atoms of a protein, [using the amino acid sequence and aligned sequences homology.](https://github.com/b0mTrady/awesome-structural-bioinformatics)

![image](https://user-images.githubusercontent.com/64801585/126504747-281b12dd-4157-4d73-a7f2-107c26494f1c.png)



* PreProcessing
  * Input Sequence 
  * Multiple Sequence Alignments
  * Structural Templates  
* Transformer (EvoFormer)
* Recycling
* Structure Module -> 3D coordinates 

![image](https://user-images.githubusercontent.com/64801585/127316142-126458b5-edf4-4bc0-8aeb-d42a24d01750.png)

![Screenshot from 2021-07-28 07-58-02](https://user-images.githubusercontent.com/64801585/127318851-d3c5f87e-75ba-4632-aa13-7b68eee2f2f8.png)

![Screenshot from 2021-07-28 07-58-54](https://user-images.githubusercontent.com/64801585/127318883-b049f5c5-9415-40b6-9de0-9eac288dcb34.png)



```python

def softmax_cross_entropy(logits, labels):
  loss = -jnp.sum(labels * jax.nn.log_softmax(logits), axis=-1)
  return jnp.asarray(loss)
  
```
If you didn't know jax's [nn.logsoftmax](https://github.com/google/jax/blob/890a41f7191fa468e2f638ba4efb9e32ad26adaa/jax/_src/nn/functions.py#L264) AF2's implemenation would not mean much to you. 

So going down the rabbit hole in Jax's nn we have the softmax function:

  (The `LogSoftmax` function, rescales elements to the range <img src="https://render.githubusercontent.com/render/math?math=(-\infty, 0)">)


```python
def log_softmax(x: Array, axis: Optional[Union[int, Tuple[int, ...]]] = -1) -> Array:  
  shifted = x - lax.stop_gradient(x.max(axis, keepdims=True))
  return shifted - jnp.log(jnp.sum(jnp.exp(shifted), axis, keepdims=True))
  ```

The accepted arguments are: 
* **x** : input array
* **axis**: the axis or axes along which the `log_softmax` should be computed. Either an integer or a tuple of integers.

and an array is returned.

Inside this function we go further down the lane to:
* [`lax.stop_gradient`](https://github.com/google/jax/blob/890a41f7191fa468e2f638ba4efb9e32ad26adaa/jax/_src/lax/lax.py#L1661) - is the identity function, that is, it returns argument `x` unchanged. However, ``stop_gradient`` prevents the flow of
  gradients during forward or reverse-mode automatic differentiation.
```python
def stop_gradient(x):
  def stop(x):
    if (dtypes.issubdtype(_dtype(x), np.floating) or
        dtypes.issubdtype(_dtype(x), np.complexfloating)):
      return ad_util.stop_gradient_p.bind(x)
    else:
      return x  # only bind primitive on inexact dtypes, to avoid some staging
  return tree_map(stop, x)
```
This in turn relies upon [`tree_map`](https://github.com/google/jax/blob/890a41f7191fa468e2f638ba4efb9e32ad26adaa/jax/_src/tree_util.py#L144)

```python 
def tree_map(f: Callable[..., Any], tree: Any, *rest: Any,
                    is_leaf: Optional[Callable[[Any], bool]] = None) -> Any:
  
  leaves, treedef = tree_flatten(tree, is_leaf)
  all_leaves = [leaves] + [treedef.flatten_up_to(r) for r in rest]
  return treedef.unflatten(f(*xs) for xs in zip(*all_leaves))
```

* `jnp.log`
* `jnp.sum`
* `jnp.exp`



[Automatic Differentiation Lecture Slides](https://www.cs.ubc.ca/~fwood/CS340/lectures/AD1.pdf)

[Gans in Jax](https://github.com/lweitkamp/GANs-JAX)

[Jax MD](https://github.com/google/jax-md)

[Open Smiles](http://opensmiles.org/opensmiles.html) - Get those old smiles out of here and imagine the wind in your hair in the drivers seat of open source smiles. The only problem - this project hasn't been updated in five years? 

