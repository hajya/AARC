AARC
====

Using residue conservation to predic protein interactions

Project proposal: In molecular biology, conservation refers to the slowness at which sections of nucleic acids or amino acids change. The more conserved something is, the more slowly it has changed throughout evolutionary history. There is a strong correlation between how (conserved) an amino acid is in a protein to the importance of the amino acid in protein-protein binding interactions. Highly conserved amino acids are more likely to occur on the surface where two proteins interact. Theoretically, by looking at evolutionary history to see how the conserved amino acids changed in one protein in relation to amino acid changes on another protein, we can glean information on how those animo acids interact.

I propose that we make a program which if given an input of two proteins, examines how conserved each amino acid is. Amino acids are scored by looking at how similar the protein sequences are in different organisms. Then the program will compare changes in amino acids between the two proteins thoughout evolutionary history. Again, by using the sequences of other current species' versions of the protein. Another score will be given, based off our measuring of how often two amino acids--one on each of the input proteins--changed concurrently. Some aspects of the project would likely include:

1) Scraping databases for protein sequences (probably using NCBI)

2) Utilizing phylogenetic trees to determine evolutionary distance between different versions of a protein

3) Developing and implementing a statistical measure of conservation (There are a lot of papers that use measures of residue conservation, so we'd have a lot of examples to use)

4) Developing and implementing a statistical measure for 'linking' two amino acids together as likely interacting. (There are some examples on how to do this already in the literature as well)


How this would be useful:

1) It would help test how conserved amino acids on interacting protein surfaces evolve concurrently in relation to another. This knowledge may be useful in protein design when engineering new proteins.

2) I work with a protease--a protein that modifies other proteins--gene of which there are more copies in fruit flies then there are in other organisms. One of these new copies is likely evolving a more specialized function in the testes, and probably modifies only a select few proteins. We have a short list of candidate proteins that it is thought to interact with, and this project would allow me to empirically test which it has more likely evolved to interact with.


Notes: 1) We would be looking at proteins sequences, like when we did the BLASTP, instead of DNA sequences. 2) I have been told by one of my professors that there is already software which determines how conserved individual amino acids are, which may be useful 3) To determine the accuracy of our program we can test pairs of proteins which are known to interact, and pairs of proteins which do not. 4) To further test the program we could compare our result to the outputs of protein-protein interaction 3D modeling software like FTDock, RDOCK, or ZDOCK. 5) Some papers which do similar, but not the same, thing: http://www.sciencemag.org/content/286/5438/295.short and http://www.proteinscience.org/details/journalArticle/112256/Physicochemical_and_residue_conservation_calculations_to_improve_the_ranking_of_.html

Also, I think this would be a novel approach. I haven't found anything yet which compares amino acid conservation on two different proteins. 
