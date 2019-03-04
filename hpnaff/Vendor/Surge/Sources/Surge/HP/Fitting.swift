//
//  Fitting.swift
//  Surge
//
//  Created by Heiko Pälike on 25/05/2015.
//  Copyright (c) 2015 Mattt Thompson. All rights reserved.
//

import Foundation
import Accelerate

/**
Calculate theta values using the normal equations.

:returns: Double array with polynomial coefficients from high to low
*/

//This matrix equation can be solved numerically, or can be inverted directly if it is well formed, to yield the solution vector
//http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
//    a=(X^(T)X)^(-1)X^(T)y.

public func polyfit(_ x: [Double], y: [Double], degree n: Int) -> [Double] {
	// 𝜃 = inverse(X' * X) * X' * y
	// Equivalent to (X' * X) * 𝜃 = X' * y hence can use la_solve

	let nCoeff = n + 1

	let m = x.count
	assert(((m - n) > 0) && nCoeff>0,"need at least y+1 points for polynomial of degree y")
	
	let yvec = la_vector_column_from_double_array(array: y)
	
	//construct Vandermonde Matrix
	var columns = [[Double]]()
	let initialColumn = [Double](repeating: 1.0 , count: m)
	
	var currentColumn = initialColumn
	
	for _ in 1...nCoeff {
		columns.append(currentColumn)
		currentColumn = mul(x, currentColumn)
	}
	
	let transpMatrix = objectFromArray(array: columns)
	let xMatrix = la_transpose(transpMatrix)
	
	let newTheta = la_solve(transpMatrix * xMatrix, transpMatrix * yvec)
	
	return newTheta.toArray().reversed()
}

public func polyfit(_ y: [Double], degree n: Int) -> [Double] {
	let x = (0..<(y.count)).map({Double($0)})
	return polyfit(x,y: y, degree: n)
}


public func detrend(_ x: [Double], y: [Double], degree n: Int = 1)->[Double]{
	let fitted = polyval(coefficients: polyfit(x,y: y,degree: n), evaluateAt: x)
	return sub(y, fitted)

}


public func detrend(_ y: [Double], degree n: Int = 1)->[Double]{
	let x = (0..<(y.count)).map({Double($0)})
	let fitted = polyval(coefficients: polyfit(y,degree: n), evaluateAt: x)
	return sub(y, fitted)
}

public func detrend(_ x: [Double], y: [Double], degree n: Int = 1)->([Double],[Double]){
	let coefficients = polyfit(x,y: y,degree: n)
	let fitted = polyval(coefficients: coefficients, evaluateAt: x)
	return (sub(y, fitted), coefficients)
}

public func detrend(_ y: [Double], degree n: Int = 1)->([Double],[Double]){
	let x = (0..<(y.count)).map({Double($0)})
	let coefficients = polyfit(y,degree: n)
	let fitted = polyval(coefficients: coefficients, evaluateAt: x)
	return (sub(y, fitted), coefficients)
}

public func detrend(real: [Double], imag: [Double], degree n: Int = 1)->((real: [Double], imag: [Double]),realCoefficients: [Double], imagCoefficients: [Double]){
    precondition(real.count == imag.count, "error in detrend of cplx: real and imag must have same number of elements.")
    //first compute magnitude and phase
    let x = (0..<(real.count)).map({Double($0)})
    let zipped = zip(real,imag)
    let magnitude = zipped.map{hypot($0.0,$0.1)}
    let phase = zipped.map{atan2($0.1,$0.0)}
    let coefficients = polyfit(magnitude,degree: n)
    let fitted = polyval(coefficients: coefficients, evaluateAt: x)
    let detrendedMagnitude = sub(magnitude, fitted)
    let zippedDetrended = zip(detrendedMagnitude,phase)
    let detrendedComplexReal = zippedDetrended.map{$0*cos($1)}
    let detrendedComplexImag = zippedDetrended.map{$0*sin($1)}
	return ((real:detrendedComplexReal,imag:detrendedComplexImag), realCoefficients: coefficients, imagCoefficients: coefficients)
}

public func detrend(_ y: [Double], coefficients: [Double])->[Double]{
    let x = (0..<(y.count)).map({Double($0)})
    let fitted = polyval(coefficients: coefficients, evaluateAt: x)
    return (sub(y, fitted))
}

/**
la_status_t to friendly string converter

:param: status The status returned from the matrix operation
*/
private func assertStatusIsSuccess(status: la_status_t) {
	switch Int32(status) {
	case LA_WARNING_POORLY_CONDITIONED:
		assert(false, "Poorly conditioned")
	case LA_INTERNAL_ERROR:
		assert(false, "Internal error")
	case LA_INVALID_PARAMETER_ERROR:
		assert(false, "Invalid parameter error")
	case LA_DIMENSION_MISMATCH_ERROR:
		assert(false, "Dimension mismatch error")
	case LA_PRECISION_MISMATCH_ERROR:
		assert(false, "Precision mismatch error")
	case LA_SINGULAR_ERROR:
		assert(false, "Singular error")
	case LA_SLICE_OUT_OF_BOUNDS_ERROR:
		assert(false, "Out of bounds error")
	case LA_SUCCESS:
		fallthrough
	default:
		break
	}
}


