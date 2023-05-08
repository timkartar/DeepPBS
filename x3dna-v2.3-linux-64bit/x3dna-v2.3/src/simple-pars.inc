    fprintf(fp, "The 'simple' base-pair and step parameters were introduced into 3DNA as of\n");
    fprintf(fp, "v2.3-2016jan01. This list of 'simple' parameters is reported by default,\n");
    fprintf(fp, "but can be turned off by specifying 'analyze -simple=false'.\n\n");

    fprintf(fp, "The simple parameters are 'intuitive' for non-Watson-Crick base pairs and\n");
    fprintf(fp, "associated base-pair steps, where the above corresponding 3DNA parameters\n");
    fprintf(fp, "may appear cryptic. Note that the following sets of simple parameters are\n");
    fprintf(fp, "for structural description only, not to be fed into the 'rebuild' program.\n");
    fprintf(fp, "Overall, they complement the rigorous characterization of base-pair\n");
    fprintf(fp, "geometry, as exemplified by the original 3DNA analyze/rebuild programs.\n\n");

    fprintf(fp, "In short, the 'simple' base-pair parameters employ the RC8--YC6 (default,\n");
    fprintf(fp, "or RN9--YN1) vector as the (long) y-axis of the pair. As before, the\n");
    fprintf(fp, "z-axis is the average of two base normals, taking consideration of the M-N\n");
    fprintf(fp, "vs M+N base-pair classification. In essence, the 'simple' parameters make\n");
    fprintf(fp, "geometrical sense through the introduction of an ad hoc base-pair\n");
    fprintf(fp, "reference frame in each case. See x3dna.org for more details. The same\n");
    fprintf(fp, "idea applies to the 'simple' inter-base-pair step parameters, which use\n");
    fprintf(fp, "consecutive C1'--C1' vectors.\n\n");

    fprintf(fp, "Here, 'angle' refers to the inter-base angle of each pair, with values in\n");
    fprintf(fp, "the range of [0, 90] degrees corresponding to the net non-planarity, i.e.,\n");
    fprintf(fp, "sqrt(buckle^2 + propeller^2).\n\n");
