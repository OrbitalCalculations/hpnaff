//
//  File.swift
//  
//
//  Created by wiggles on 06/12/2021.
//

import Foundation
import DoggieMath

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
  assert(((m - n) > 0) && (nCoeff > 0),"need at least y+1 points for polynomial of degree y")
  
  let poly = MyPolynomial.fit(x: x, y: y, degree: n)
  return poly.coefficients//.reversed()
}

public func polyfit(_ y: [Double], degree n: Int) -> [Double] {
  let x = (0..<(y.count)).map({Double($0)})
  return polyfit(x,y: y, degree: n)
}

// MARK: Polynomials

public func polyval(coefficients: [Double], evaluateAt x: [Double])->[Double]{
  //var results = [Double](repeating: 0.0, count: x.count)
  let polynomialDegree = coefficients.count - 1
  assert(polynomialDegree >= 0)
  let poly = MyPolynomial(coefficients)
  let results = x.map{poly.evaluated($0)}
  //poly.evaluated(<#T##Double#>)
  //vDSP_vpolyD(coefficients, 1, x, 1, &results, 1, vDSP_Length(x.count), vDSP_Length(polynomialDegree))
  return results
}


public func detrend(_ x: [Double], y: [Double], degree n: Int = 1)->[Double]{
  let fitted = polyval(coefficients: polyfit(x, y: y,degree: n), evaluateAt: x)
  return zip(y,fitted).map{$0-$1}

}


public func detrend(_ y: [Double], degree n: Int = 1)->([Double],[Double]){
  let x = (0..<(y.count)).map({Double($0)})
  let coefficients = polyfit(y,degree: n)
  let fitted = polyval(coefficients: coefficients, evaluateAt: x)
  return (zip(y,fitted).map{$0-$1}, coefficients)
}

public func detrend(_ x: [Double], y: [Double], degree n: Int = 1)->([Double],[Double]){
  let coefficients = polyfit(x,y: y,degree: n)
  let fitted = polyval(coefficients: coefficients, evaluateAt: x)
  return (zip(y,fitted).map{$0-$1}, coefficients)
}

//public func detrend(_ y: [Double], degree n: Int = 1)->([Double],[Double]){
//  let x = (0..<(y.count)).map({Double($0)})
//  let coefficients = polyfit(y,degree: n)
//  let fitted = polyval(coefficients: coefficients, evaluateAt: x)
//  return (zip(y,fitted).map{$0-$1}, coefficients)
//}

public func detrend(real: [Double], imag: [Double], degree n: Int = 1)->((real: [Double], imag: [Double]),realCoefficients: [Double], imagCoefficients: [Double]){
  precondition(real.count == imag.count, "error in detrend of cplx: real and imag must have same number of elements.")
  //first compute magnitude and phase
  let x = (0..<(real.count)).map({Double($0)})
  let zipped = zip(real,imag)
  let magnitude = zipped.map{hypot($0.0,$0.1)}
  let phase = zipped.map{atan2($0.1,$0.0)}
  let coefficients = polyfit(x, y: magnitude,degree: n)
  let fitted = polyval(coefficients: coefficients, evaluateAt: x)
  let detrendedMagnitude = zip(magnitude,fitted).map{$0-$1}
  let zippedDetrended = zip(detrendedMagnitude,phase)
  let detrendedComplexReal = zippedDetrended.map{$0*cos($1)}
  let detrendedComplexImag = zippedDetrended.map{$0*sin($1)}
  return ((real:detrendedComplexReal,imag:detrendedComplexImag), realCoefficients: coefficients, imagCoefficients: coefficients)
}

public func detrend(_ y: [Double], coefficients: [Double])->[Double]{
  let x = (0..<(y.count)).map({Double($0)})
  let fitted = polyval(coefficients: coefficients, evaluateAt: x)
  return zip(y,fitted).map{$0-$1}
}
