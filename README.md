**Inspiration**
Antibiotic resistance is increasingly dangerous to public health and is one of the major concerns in drug discovery. It is therefore very important to find new drugs in a highly time and cost-effective manner to keep up with the rate of the antibiotic resistance. Virtually synthesizing potentially active compounds is quite exciting, and crucial in the drug discovery research. 

**What it does**
The program takes a chemical compound and generates analogs of that target compound. It does so by fragmenting into synthetically accessible building blocks and replacing each with other building blocks with similar chemical properties. Essentially, it generates a virtual chemical library by assembling chemically similar fragments and generating analogs of the target compound.

**How I built it**
It was all built in python, using RDKit library.
First, I built a chemical database which contains experimentally synthesizable 10k molecules from ZINC15 database (https://zinc15.docking.org/). Then, all the chemicals from that database are fragmented into synthetic building blocks to make a ‘chemical fragments’ database. It hosts unique substructures derived from ZINC15 chemicals and serves as a database for substructure replacement.
As the program input, ChemX takes a smile string (digital format of a chemical compound) and converts it into a mol object. Then, it performs step-by-step chemical manipulations to fragment the compound into synthetic building blocks. For each building block, the program searches for fragments in the ‘chemical fragments’ database which share similar chemical properties (about ~200 2d descriptors can be computed for each chemical). Similarity was done by computing the Euclidean distance between lists containing chemical properties of compounds. The fragments from “chemical fragments” database with the highest similarity scores for each fragment identified in the target drug are gathered. Finally, ChemX generates different arrangements of those chosen fragments to build new virtual compounds.

**Challenges I ran into**
Researching the right database and tools to use and manipulating the backends of different libraries to fit the needs of the program were challenging. Building fragment database, parsing and joining them to make valid molecules is also stressful. Handling a lot of chemistry toolkit related issues such as invalid valencies, bond counts, atoms and standardizing the structures was an issue as well.

**Accomplishments that I'm proud of**
I successfully prepared a prototype for the first set of functionalities that I aimed for ChemX, which was to generate virtual chemical libraries of analogs by applying the fragment-replacement approach. I generated sample chemical analog libraries for antimalarial drugs including mefloquine and quinine.
I am especially proud that I have a working prototype at the end of the day, especially when I don’t have any teammate. 

** What I learned**
I learned a lot about chemical manipulations such as fragmenting, merging and substituting substructures. I also read a lot about molecular graphs and applications of generative models to create molecules, though I didn’t end up implementing them as part of this project.

** What's next for ChemX**
I plan to research and implement better approaches for merging the fragments together, explore different methods to find similar substructures and compare results and optimize the generated chemical library by incorporating more predictors/features such as toxicity, synthetic feasibilities. There is certainly additional work to do regarding the curation of the compounds, removing duplicates, joining fragments in a more controlled and restricted manner, etc. 
