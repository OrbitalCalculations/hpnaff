//
//  LinearRegression.swift
//  LinearAlgebraExtensions
//
//  Created by Damien Pontifex on 16/09/2014.
//  Copyright (c) 2014 Damien Pontifex. All rights reserved.
//

//        var x : [Double] = [0.0, 1.0, 2.0]
//let y : [Double] = [2.0, 3.0, 6.0]


//let lr = LinearRegression(x,y)
//let res = lr.normalEquations()

//println(res)


import Foundation
import Accelerate

/**
*  Class to calculate linear regression of a set of data points
*/
public class LinearRegression {
	
	private var m: Int
	private var x: la_object_t
	private var y: la_object_t
	private var theta: la_object_t
	private var alpha: Double
	private var numIterations: Int
	
	public init(_ x: la_object_t, _ y: la_object_t, theta: la_object_t? = nil, alpha: Double = 0.01, numIterations: Int = 1500) {
		self.m = Int(la_matrix_rows(x))
		//let cols = Int(la_matrix_cols(x))
		
		let ones = la_ones_matrix(rows: Int(la_matrix_rows(x)))
		self.x = x.prependColumnsFrom(ones)
		self.y = y
		
		if theta != nil {
			self.theta = theta!
		} else {
			self.theta = la_zero_matrix(columns: Int(la_matrix_cols(self.x)))
		}
		self.alpha = alpha
		self.numIterations = numIterations
	}
	
	/**
	Initialise the LinearRegression object
	
	:param: x             Input data values
	:param: y             Expected output values
	:param: theta         Initial theta value
	:param: alpha         Learning rate
	:param: numIterations Number of iterations to run the algorithm
	
	:returns: Instance of the LinearRegression
	*/
	public init(_ x: [Double], _ y: [Double], theta: [Double]? = nil, alpha: Double = 0.01, numIterations: Int = 1500) {
		self.m = x.count
		
		var xValues = Array<Double>(repeating: 1.0, count: m * 3)
		for i in stride(from: 0, to: x.count * 3, by: 3)  { //stride(from: 0, to: x.count * 3, by: 3) {
			xValues[i+1] = x[i / 3]
		}
        for i in stride(from: 0, to: x.count * 3, by: 3)  { //stride(from: 0, to: x.count * 3, by: 3) {
            xValues[i+2] = x[i / 3] * x[i / 3]
        }
		
		let xMatrix = la_matrix_from_double_array(xValues, rows: m, columns: 3)
		
		self.x = xMatrix
		
		self.y = la_vector_column_from_double_array(array: y)
		var thetaArray = theta
		if thetaArray == nil {
			// Default the theta array to zeros with appropriate size of x array
			thetaArray = Array<Double>(repeating: 0.0, count: Int(3))
		}
		self.theta = la_vector_row_from_double_array(array: thetaArray!)
		self.alpha = alpha
		self.numIterations = numIterations
	}
	
	/**
	Run the gradient descent algorithm for linear regression
	
	:param: returnCostHistory Boolean value to indicate whether the cost history should be returned. Some performance hit if this is true due to doing array computations in every iteration rather than just at the end
	
	:returns: Tuple of theta values and optional jHistory parameter if requested
	*/
	public func gradientDescent(returnCostHistory: Bool = false) -> (theta: [Double], jHistory: [Double]?) {
		
		print("initial theta \(theta.description)")
		print(x.description)
		print(y.description)
		
		// Number of training examples
		let alphaOverM = alpha / Double(m)
		
		var jHistory: [Double]?
		if returnCostHistory {
			jHistory = Array<Double>(repeating: 0.0, count: numIterations)
		}
		
		for iter in 0..<numIterations {
			// h(x) = transpose(theta) * x
			let prediction = x * theta
			let errors = prediction - y
			let sum = la_transpose(errors) * x
			let partial = la_transpose(sum) * alphaOverM
			
			// Simultaneous theta update:
			// theta_j = theta_j - alpha / m * sum_{i=1}^m (h_theta(x^(i)) - y^(i)) * x_j^(i)
			theta = theta - partial
			
			print("next theta \(theta.description)")
			
			if returnCostHistory {
				jHistory![iter] = computeCost()
			}
		}
		
		let thetaArray = theta.toArray()
		
		return (theta: thetaArray, jHistory: jHistory)
	}
	
	/**
	Calculate theta values using the normal equations.
	
	:returns: Double array with theta coefficients
	*/
    
    //This matrix equation can be solved numerically, or can be inverted directly if it is well formed, to yield the solution vector
    //http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
//    a=(X^(T)X)^(-1)X^(T)y.
	public func normalEquations() -> [Double] {
		// 𝜃 = inverse(X' * X) * X' * y
		// Equivalent to (X' * X) * 𝜃 = X' * y hence can use la_solve
		let newTheta = la_solve(la_transpose(x) * x, la_transpose(x) * y)
		
		return newTheta.toArray()
	}
	
	/**
	Computes the cost with our current values of theta
	
	:returns: Cost value for the data with value of theta
	*/
	public func computeCost() -> Double {
		
		let twoM = 2.0 * Double(m)
		
		let xTheta = x * theta
		let diff = xTheta - y
		
		var partial = diff^2
		partial = partial / twoM
		
		return partial.toArray().first!
	}
}
