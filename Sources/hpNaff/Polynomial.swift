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

/// Represents a polynomial function, e.g. `2 + 3x + 4x²`.
public struct MyPolynomial: Codable, Equatable {
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

//  public func callAsFunction(_ temperature: Temperature) -> Double {
//    evaluated(temperature.kelvin)
//  }

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

extension MyPolynomial: ExpressibleByArrayLiteral {
  public init(arrayLiteral elements: Double...) {
    self.coefficients = elements
  }
}

extension MyPolynomial: CustomStringConvertible {
  public var description: String {
    var s: String = ""
    for (i, c) in coefficients.enumerated() {
      s += "c\(i): \(String(format: "%.6e", c))"
    }
    return s
  }
}

extension MyPolynomial {
  public static func fit(x: [Double], y: [Double], degree n: Int = 5) -> MyPolynomial {
    /// degree of polynomial to fit the data
    var n: Int = n
    /// no. of data points
    let N: Int = min(x.count, y.count)

    var X: [Double] = Array(repeating: 0.0, count: 2 * n + 1)

    for i in X.indices {
      for j in 0..<N {
        // consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
        X[i] += pow(x[j], Double(i))
      }
    }
    var a: [Double] = Array(repeating: 0, count: n + 1)
    /// B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
    var B: [[Double]] = Array(repeating: Array(repeating: 0, count: n + 2), count: n + 1)

    for i in 0...n {
      for j in 0...n {
        // Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
        B[i][j] = X[i + j]
      }
    }

    /// Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
    var Y: [Double] = Array(repeating: 0, count: n + 1)

    for i in 0..<(n + 1) {
      Y[i] = 0
      for j in 0..<N {
        // consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
        Y[i] += pow(x[j], Double(i)) * y[j]
      }
    }

    for i in 0...n {
      // load the values of Y as the last column of B(Normal Matrix but augmented)
      B[i][n + 1] = Y[i]
    }

    n += 1
    for i in 0..<n {
      // From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
      for k in (i + 1)..<n {
        if B[i][i] < B[k][i] {
          for j in 0...n {
            let temp = B[i][j]
            B[i][j] = B[k][j]
            B[k][j] = temp
          }
        }
      }
    }

    for i in 0..<(n - 1) {  // loop to perform the gauss elimination
      for k in (i + 1)..<n {
        let t = B[k][i] / B[i][i]
        for j in 0...n {
          // make the elements below the pivot elements equal to zero or elimnate the variables
          B[k][j] -= t * B[i][j]
        }
      }
    }

    for i in (0..<(n - 1)).reversed() { // back-substitution
      // x is an array whose values correspond to the values of x,y,z..
      // make the variable to be calculated equal to the rhs of the last equation
      a[i] = B[i][n]
      for j in 0..<n {
        if j != i {
          // then subtract all the lhs values except the coefficient of the variable whose value is being calculated
          a[i] -= B[i][j] * a[j]
        }
      }
      a[i] /= B[i][i] // now finally divide the rhs by the coefficient of the variable to be calculated
    }
    a.removeLast()
    return self.init(a)
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
