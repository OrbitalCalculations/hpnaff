//
//  Power2.swift
//  Surge
//
//  Created by Heiko Pälike on 26/05/2015.
//  Copyright (c) 2015 Mattt Thompson. All rights reserved.
//

import Accelerate


// MARK: Polynomials
public func polyval(coefficients: [Float], evaluateAt x: [Float])->[Float]{
	var results = [Float](repeating: 0.0, count: x.count)
	let polynomialDegree = coefficients.count - 1
	assert(polynomialDegree >= 0)
	vDSP_vpoly(coefficients, 1, x, 1, &results, 1, vDSP_Length(x.count), vDSP_Length(polynomialDegree))
	return results
}

public func polyval(coefficients: [Double], evaluateAt x: [Double])->[Double]{
	var results = [Double](repeating: 0.0, count: x.count)
	let polynomialDegree = coefficients.count - 1
	assert(polynomialDegree >= 0)
	vDSP_vpolyD(coefficients, 1, x, 1, &results, 1, vDSP_Length(x.count), vDSP_Length(polynomialDegree))
	return results
}
