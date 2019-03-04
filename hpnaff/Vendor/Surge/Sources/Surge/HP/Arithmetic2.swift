// Arithmetic.swift
//
// Copyright (c) 2014–2015 Mattt Thompson (http://mattt.me)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
import Foundation
import Accelerate
// MARK: Complex Vector Amplitude

public func amp(real:[Float], imag:[Float])->[Float]{
	var result = [Float](repeating: 0.0, count: real.count)
	
	var real = real
	var imag = imag
	
	var splitcomplex = DSPSplitComplex(realp: &real, imagp: &imag)
	vDSP_zvabs(&splitcomplex, 1, &result, 1, vDSP_Length(result.count))
	return result
}

public func amp(real:[Double], imag:[Double])->[Double]{
	var result = [Double](repeating: 0.0, count: real.count)
	
	var real = real
	var imag = imag
	
	var splitcomplex = DSPDoubleSplitComplex(realp: &real, imagp: &imag)
	vDSP_zvabsD(&splitcomplex, 1, &result, 1, vDSP_Length(result.count))
	return result
}

// MARK: Subtract
/*
public func sub(_ x: [Float], _ y: [Float]) -> [Float] {
	var results = [Float](repeating: 0.0, count: x.count)
	vDSP_vsub(y, 1, x, 1, &results, 1, vDSP_Length(x.count))
	return results
}

public func sub(_ x: [Double], _ y: [Double]) -> [Double] {
	var results = [Double](repeating: 0.0, count: x.count)
	vDSP_vsubD(y, 1, x, 1, &results, 1, vDSP_Length(x.count))
	return results
}*/


// MARK: Vector Sub Multiply
// (a-b)*z

public func submul(_ x: [Float], _ y: [Float], z: [Float]) -> [Float] {
	var results = [Float](repeating: 0.0, count: x.count)
	vDSP_vsbm(x, 1, y, 1, z, 1, &results, 1, vDSP_Length(x.count))
	return results
}


public func submul(_ x: [Double], _ y: [Double], _ z: [Double]) -> [Double] {
	var results = [Double](repeating: 0.0, count: x.count)
	vDSP_vsbmD(x, 1, y, 1, z, 1, &results, 1, vDSP_Length(x.count))
	return results
}

// MARK: Vector Scalar Add

public func vsadd(_ x: [Float], _ c: Float)->[Float]{
	var results = [Float](repeating: 0.0, count: x.count)
	vDSP_vsadd(x, 1, [c], &results, 1, vDSP_Length(x.count))
	return results
}

public func vsadd(_ x: [Double], _ c: Double)->[Double]{
	var results = [Double](repeating: 0.0, count: x.count)
	vDSP_vsaddD(x, 1, [c], &results, 1, vDSP_Length(x.count))
	return results
}

// MARK: Vector Scalar Multiply
public func vsmul(_ x: [Float], _ c: Float)->[Float]{
	var results = [Float](repeating: 0.0, count: x.count)
	vDSP_vsmul(x, 1, [c], &results, 1, vDSP_Length(x.count))
	return results
}

public func vsmul(_ x: [Double], _ c: Double)->[Double]{
	var results = [Double](repeating: 0.0, count: x.count)
	vDSP_vsmulD(x, 1, [c], &results, 1, vDSP_Length(x.count))
	return results
}
// MARK: Vector squared
public func sqr(_ x: [Float])->[Float]{
    var results = [Float](repeating: 0.0, count: x.count)
    vDSP_vsq(x, 1, &results, 1, vDSP_Length(x.count))
    return results
}

public func sqr(_ x: [Double])->[Double]{
    var results = [Double](repeating: 0.0, count: x.count)
    vDSP_vsqD(x, 1, &results, 1, vDSP_Length(x.count))
    return results
}

// MARK: Vector Scalar squared sum
public func svesq(_ x: [Float])->Float{
	var result : Float = 0.0
	vDSP_svesq(x, 1, &result, vDSP_Length(x.count))
	return result
}

public func svesq(_ x: [Double])->Double{
	var result : Double = 0.0
	vDSP_svesqD(x, 1, &result, vDSP_Length(x.count))
	return result
}
