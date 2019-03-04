//
//  Statistics.swift
//  Surge
//
//  Created by Heiko Pälike on 02/06/2015.
//  Copyright (c) 2015 Mattt Thompson. All rights reserved.
//

import Foundation
//import Surge

/*
rnorm(npts,sd=sqrt(noisevar))

{"rnorm",	do_random2,	8,	11,	3,	{PP_FUNCALL, PREC_FN,	0}},
*/
/**
	Mersenne Twister random number generator used to compute distribution

	:param: **n**         Number of values to generate and return.
	:param: **mean**      Desired Mean of results (default = 0.0).
	:param: **sd**        Desired Standard deviation (default = 1.0, calculated with n-1)

	:returns: Random distribution with desired mean and standard deviation.

*/
public func rnorm(n: Int, rnormmean: Double = 0, sd: Double = 1.0)->[Double]{

    var randoms = [Double]()
    
    let mtrandom = MTRandom()!
    
    for _ in 0..<n {
        randoms.append(mtrandom.randomDouble())
    }
    
    let initial_mean = mean(randoms)  //Surge.mean(randoms)
    randoms = vsadd(randoms, -1.0 * initial_mean)
    let initial_sd = sqrt(svesq(randoms)/(Double(n) - 1.0))
    
    let normalized_randoms = vsmul(randoms, (1.0 / initial_sd))
    
    let results = vsadd(vsmul(normalized_randoms, sd), rnormmean)
    
    return results
}
