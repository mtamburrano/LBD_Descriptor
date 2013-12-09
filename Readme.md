The source code of the proposed line matching algrothim is  included. The list of files is:

LineStructure.hh: the definition of the line structure;
EDLineDetector.cpp and EDLineDetector.hh: the class of EDLineDetector which is our implementation of EDLine detector;
LineDescriptor.cpp and LineDescriptor.cpp: the class of LBD line descriptor which includes extract lines in scaple space and construct their LBD descriptors;
PairwiseLineMatching.cpp and PairwiseLineMatching.hh: the class of line matching based on the local appearance similarities and pairwise geometric consistencies.
TestLineMatchingAlgorithm.cpp: an example to use the above classes.

CMakeLists.txt: the cmake file to compile the project.

references:

@article{Zhang2013794,
title = "An efficient and robust line segment matching approach based on \{LBD\} descriptor and pairwise geometric consistency ",
journal = "Journal of Visual Communication and Image Representation ",
year = "2013",
issn = "1047-3203",
doi = "http://dx.doi.org/10.1016/j.jvcir.2013.05.006",
url = "http://www.sciencedirect.com/science/article/pii/S1047320313000874",
author = "Lilian Zhang and Reinhard Koch"
}

