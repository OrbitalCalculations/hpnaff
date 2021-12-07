//
//  Copyright 2017 Daniel Müllenborn
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// from: https://github.com/damuellen/Utilities


import Foundation
import Numerics

public enum PolyfitError: String, Error {
  case SingularMatrix
  case TooFewDatapointsForOrder
  case OrderGreaterThanMaxOrderFive
}

/// Represents a polynomial function, e.g. `2 + 3x + 4x²`.
public struct Poly: Codable, Equatable {
  /// Represents the coefficients of the polynomial
  public let coefficients: [Double]

  public init(coeffs: Double...) {
    self.coefficients = coeffs
  }

  public init(_ array: [Double]) {
    self.coefficients = array
  }

  public var indices: CountableRange<Int> { coefficients.indices }

  public var isEmpty: Bool { coefficients.isEmpty }

  public var isInapplicable: Bool { coefficients.count < 2 }

  @_transparent func evaluated(_ value: Double) -> Double {
    // Use Horner’s Method for solving
    coefficients.reversed().reduce(into: 0.0) { result, coefficient in
      result = coefficient.addingProduct(result, value)
    }
  }

  public func callAsFunction(_ value: Double) -> Double {
    evaluated(value)
  }

  public func callAsFunction(_ ratio: Ratio) -> Double {
    evaluated(ratio.quotient)
  }

  public subscript(index: Int) -> Double {
    coefficients[index]
  }
}

extension Poly: ExpressibleByArrayLiteral {
  public init(arrayLiteral elements: Double...) {
    self.coefficients = elements
  }
}

extension Poly: CustomStringConvertible {
  public var description: String {
    var s: String = ""
    for (i, c) in coefficients.enumerated() {
      s += "c\(i): \(String(format: "%.6e", c))"
    }
    return s
  }
}

extension Poly {
  //https://github.com/natedomin/polyfit
  //----------------------------------------------------
  //
  // METHOD:  polyfit
  //
  // INPUTS:  dependentValues[0..(countOfElements-1)]
  //          independentValues[0...(countOfElements-1)]
  //          countOfElements
  //          order - Order of the polynomial fitting
  //
  // OUTPUTS: coefficients[0..order] - indexed by term
  //               (the (coef*x^3) is coefficients[3])
  //
  //----------------------------------------------------
  //int polyfit(const double* const dependentValues,
  //            const double* const independentValues,
  //            unsigned int        countOfElements,
  //            unsigned int        order,
  //            double*             coefficients)
  static func fit(x dependentValues: [Double], y independentValues: [Double], order: Int) throws -> Poly {
    // Declarations...
    // ----------------------------------
    let maxOrder = 5
    var B = [Double](repeating: 0.0, count: order + 1)
    var P = [Double](repeating: 0.0, count: ((order+1) * 2)+1)
    var A = [Double](repeating: 0.0, count: (order + 1)*2*(order + 1))
    var coefficients = [Double](repeating: 0.0, count: order + 1)

    // Verify initial conditions....
    // ----------------------------------

    // This method requires that the countOfElements >
    // (order+1)
                                          
    let countOfElements = dependentValues.count
    guard (countOfElements > order) else {
      throw PolyfitError.TooFewDatapointsForOrder
    }

    // This method has imposed an arbitrary bound of
    // order <= maxOrder.  Increase maxOrder if necessary.
    guard (order <= maxOrder) else {
      throw PolyfitError.OrderGreaterThanMaxOrderFive
    }
    
    // Begin Code...
    // ----------------------------------

    // Identify the column vector
    for ii in 0..<countOfElements {
      let x    = dependentValues[ii]
      let y    = independentValues[ii]
      var powx = 1.0

      for jj in 0..<(order + 1) {
        B[jj] = B[jj] + (y * powx)
        powx *= x
      }
    }
    // Initialize the PowX array
    P[0] = Double(countOfElements)

    // Compute the sum of the Powers of X
    for ii in 0..<countOfElements {
      let x    = dependentValues[ii]
      var powx = dependentValues[ii]

      for jj in 1 ..< ((2 * (order + 1)) + 1) {
            P[jj] = P[jj] + powx;
            powx  *= x
        }
    }

    // Initialize the reduction matrix
    //
    for ii in 0..<(order + 1) {
      for jj in 0..<(order + 1) {
        A[(ii * (2 * (order + 1))) + jj] = P[ii+jj];
      }
      A[(ii*(2 * (order + 1))) + (ii + (order + 1))] = 1.0
    }

    // Move the Identity matrix portion of the redux matrix
    // to the left side (find the inverse of the left side
    // of the redux matrix
    for ii in 0..<(order + 1) {
      let x = A[(ii * (2 * (order + 1))) + ii]
      if (x != 0.0) {
        for kk in 0..<(2 * (order + 1)) {
          A[(ii * (2 * (order + 1))) + kk] =
                    A[(ii * (2 * (order + 1))) + kk] / x
        }
        for jj  in 0..<(order + 1) {
          if ((jj - ii) != 0) {
            let y = A[(jj * (2 * (order + 1))) + ii]
            for kk in 0..<(2 * (order + 1)) {
              A[(jj * (2 * (order + 1))) + kk] =
                  A[(jj * (2 * (order + 1))) + kk] -
                  y * A[(ii * (2 * (order + 1))) + kk]
            }
          }
        }
      } else {
      // Cannot work with singular matrices
        throw PolyfitError.SingularMatrix
      }
    }

    // Calculate and Identify the coefficients
    for ii in 0..<(order + 1) {
      var x = 0.0
      for kk in 0..<(order + 1) {
        x += (A[(ii * (2 * (order + 1))) + (kk + (order + 1))] * B[kk])
      }
      coefficients[ii] = x;
    }
    return self.init(coefficients)
  }
}

public struct Ratio: CustomStringConvertible, Codable {
  public var quotient: Double

  public var isZero: Bool { self == .zero }

  public static var zero: Ratio { Ratio(0) }

  public var percentage: Double { quotient * 100.0 }

  public var description: String {
    String(format: "%3.1f", percentage) + "%"
  }

  public init(percent: Double) {
    self.quotient = percent / 100
  }

  public init(_ value: Double) {
    precondition(0...1.01 ~= value, "Ratio out of range.")
    self.quotient = value > 1 ? 1 : value
  }

  public init(_ value: Double, cap: Double) {
    precondition(0 <= value, "Ratio out of range.")
    self.quotient = min(value, cap)
  }

  public mutating func limited(to max: Ratio) {
    quotient = min(max.quotient, quotient)
  }
}

extension Ratio: ExpressibleByFloatLiteral {
  public init(floatLiteral value: Double) {
    self.quotient = value
  }
}

extension Ratio: Equatable {
  public static func == (lhs: Ratio, rhs: Ratio) -> Bool {
    lhs.quotient == rhs.quotient
  }
}

extension Ratio: Comparable {
  public static func < (lhs: Ratio, rhs: Ratio) -> Bool {
    lhs.quotient < rhs.quotient
  }
}

extension Ratio {
  public var multiBar: String {
    let (bar_chunks, remainder) = Int(quotient * 80)
      .quotientAndRemainder(dividingBy: 8)
    let full = UnicodeScalar("█").value
    let fractionalPart = remainder > 0
      ? String(UnicodeScalar(full + UInt32(8 - remainder))!) : ""
    return String(repeating: "█", count: bar_chunks)
      + fractionalPart
      + String(repeating: " ", count: 10 - bar_chunks)
      + description
  }

  public var singleBar: String {
    let bar = Int(quotient * 7)
    let full = UnicodeScalar("█").value
    let block = String(UnicodeScalar(full - UInt32(7 - bar))!)
    return block + " " + description
  }
}

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
  
  do {
    return try Poly.fit(x: x, y: y, order: n).coefficients
  } catch {
    fatalError("Error in obtaining coefficients from Polyfit: \(error)")
  }
}

public func polyfit(_ y: [Double], degree n: Int) -> [Double] {
  let x = (0..<(y.count)).map({Double($0)})
  return polyfit(x,y: y, degree: n)
}

// MARK: Polynomials

public func polyval(coefficients: [Double], evaluateAt x: [Double]) -> [Double] {
  let polynomialDegree = coefficients.count - 1
  assert(polynomialDegree >= 0)
  let poly = Poly(coefficients)
  return x.map{poly.evaluated($0)}
}


public func detrend(_ x: [Double], y: [Double], degree n: Int = 1)->[Double]{
  let fitted = polyval(coefficients: polyfit(x, y: y,degree: n), evaluateAt: x)
  return zip(y,fitted).map{$0-$1}

}


public func detrend(_ y: [Double], degree n: Int = 1) -> ([Double],[Double]) {
  let x = (0..<(y.count)).map({Double($0)})
  let coefficients = polyfit(y,degree: n)
  let fitted = polyval(coefficients: coefficients, evaluateAt: x)
  return (zip(y,fitted).map{$0-$1}, coefficients)
}

public func detrend(_ x: [Double], y: [Double], degree n: Int = 1) -> ([Double],[Double]) {
  let coefficients = polyfit(x,y: y,degree: n)
  let fitted = polyval(coefficients: coefficients, evaluateAt: x)
  return (zip(y,fitted).map{$0-$1}, coefficients)
}


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
