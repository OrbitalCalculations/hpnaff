//
//  naff_math.swift
//  hpnaff
//
//  Created by Heiko Pälike on 29/04/2015.
//  Copyright (c) 2015 Heiko Pälike. All rights reserved.
//

import Foundation
import DoggieMath
import Numerics
//#if canImport(Accelerate)
//	import Accelerate
//#endif
//import Surge


extension Array where Element: Real {
  @inlinable
  public func squared() -> Self {
    return self.map { $0*$0 }
  }
}

extension Array where Element: Real {
  @inlinable
  public func sum() -> Self.Element {
    return self.reduce(0, +)
  }
}


extension Real {
  @inlinable
  public func svesq(_ arg: [Double]) -> Double {
    return arg.reduce(0.0, {$0 + $1*$1})
  }
}

@inlinable
public func sincos(_ arg: [Double]) -> ([Double], [Double]) {
  let sine = arg.map(sin)
  let cosine = arg.map(cos)
  return (sine, cosine)
}

func NAFFFunc(NAFFData : [Double], NAFFdt: Double, omega: Double) -> Double {
	let range = 0 ..< NAFFData.count
  let omegaterm = range.map { Double($0)*omega*NAFFdt }
	let (sine, cosine) = sincos(omegaterm)
  let sum1 = zip(cosine, NAFFData).reduce(0, {$0 + $1.0*$1.1})
  let sum2 = zip(sine, NAFFData).reduce(0, {$0 + $1.0*$1.1})
  return sum1 * sum1 + sum2 * sum2
}


func NAFFFunc2(NAFFData : [Double], NAFFCplxData : [Double], NAFFdt: Double, omega: Double) -> Double {
	let range = 0..<NAFFData.count
  let omegaterm = range.map { Double($0)*omega*NAFFdt }
  let (sine, cosine) = sincos(omegaterm)
  var sum1: Double = zip(cosine, NAFFData).reduce(0, {$0 + $1.0*$1.1})
  sum1 += zip(sine, NAFFCplxData).reduce(0, {$0 + $1.0*$1.1}) //real
  var sum2: Double = zip(cosine, NAFFCplxData).reduce(0, {$0 + $1.0*$1.1})
  sum2 -= zip(sine, NAFFData).reduce(0, {$0 + $1.0*$1.1})  //imag
  return sum1 * sum1 + sum2 * sum2
}

public func sqr(x: Int)->Int{
	return x*x
}

public func hanning(points: Int)->[Double]{
	assert(points>1, "Error in Hanning Function: Number of points should be > 1")
	//#if canImport(Accelerate)
	//var result = [Double](repeating: 0.0, count: points)
	//vDSP_hann_windowD(&result,vDSP_Length(points),Int32(vDSP_HANN_DENORM))
	//return result
	//#else
  return (0..<points).lazy.map{ (1.0 - cos(2.0*Double.pi*Double($0)/(Double(points))))/2.0 }
	//#endif
}

public struct Hanning{
	private static var _points : Int = 0
	private static var _data = [Double]()
	
	public static var points : Int {get{return Hanning._points}}
	public static var data : [Double] {get{return Hanning._data}}
	
	init(points:Int){
		if (points != Hanning.points){
			Hanning._data = hanning(points: points)
			Hanning._points = points
		}
	}
}

public func adjustFrequencyHalfPlane(frequency: Double, phase0: Double, phase1: Double, dt: Double)->Double{
	var phase0 : Double = phase0
	var phase1 : Double = phase1
	if (fabs(phase0-phase1)>Double.pi) {
		if (phase0 < phase1){
			phase0 += Double.pi*2.0
		}
		else{
			phase1 += Double.pi*2.0
		}
	}
	if (phase0 < phase1){
		return frequency
	}
	return 1.0/dt - frequency;
}

func CalculatePhaseAndAmplitudeFromFreq(
  hanning: [Double],
  data: [Double],
  nPoints: Int,
  dt: Double,
  frequency: Double,
  t0: Double) ->
    (newData:[Double], amplitude:Double, phase:Double, significance:Double) {
	
	assert (nPoints>1,"number of Points in CalculatePhaseAndAmplitudeFromFreq must be >1")
	assert (nPoints == hanning.count,"number of hanning window points must equal data points")
	
	let freq0 : Double = frequency  // * 2.0 * M_PI
	
	let range = 0 ..< nPoints
  let sincosarg = range.map {Double($0) * freq0 * dt}
	let (sine, cosine) = sincos(sincosarg)

  let sum_ee1 = zip(cosine, hanning).reduce(0.0, {$0 + $1.0*$1.0*$1.1})
  let sum_ee2 = zip(  sine, hanning).reduce(0.0, {$0 + $1.0*$1.0*$1.1})
  
  /* these are the overlap sums */
  let sum_ef1 = zip(cosine, data).reduce(0.0, {$0 + $1.0*$1.1})
  let sum_ef2 = zip(  sine, data).reduce(0.0, {$0 + $1.0*$1.1})
	
  let sum1 = data.reduce(0.0, {$0 + $1*$1})// svesq(data) //Accelerate/Surge
    //    NAFFData[i] -= (sum_ef1/sum_ee1*cosine[i] + sum_ef2/sum_ee2*sine[i])*hanning[i];

  let term1 = zip(cosine, hanning).lazy.map {$0 * sum_ef1/sum_ee1 * $1}
  let term2 = zip(  sine, hanning).lazy.map {$0 * sum_ef2/sum_ee2 * $1 }
  let term = zip(term1, term2).lazy.map{$0 + $1}
  let newData = zip(data, term).map{$0 - $1} //sub(data, term) //Accelerate/Surge
	
	let sum2 = newData.reduce(0.0, {$0 + $1*$1})
	
	let significance = (sum1>0.0) ? (sum2/sum1) : (-1.0)
	
	let freq0_2 = frequency / (Double.pi*2.0)
	
	let amplitude : Double = sqrt(pow(sum_ef1/sum_ee1, 2.0) + pow(sum_ef2/sum_ee2, 2.0))
	/* compute the phase and ensure that it is in the range [-PI, PI] */
	var phase = fmod(atan2(-sum_ef2/sum_ee2, sum_ef1/sum_ee1) +
		freq0_2*t0*Double.pi*2.0, Double.pi * 2.0);
	
	if (phase < -Double.pi){
		phase += Double.pi * 2.0
	}
	if (phase > Double.pi){
		phase -= Double.pi * 2.0
	}
	return (newData,amplitude,phase,significance)

}
func CalculatePhaseAndAmplitudeFromFreqCplx2(hanning: [Double], real: [Double], imag: [Double], nPoints: Int, dt: Double, frequency: Double, t0: Double)->(newReal:[Double], newImag: [Double], amplitude:Double, phase:Double,significance:Double){
	
	assert (nPoints>1,"number of Points in CalculatePhaseAndAmplitudeFromFreq must be >1")
	assert (nPoints == hanning.count,"number of hanning window points must equal data points")
	
	let freq0 : Double = frequency  // * 2.0 * M_PI
	
	let range = 0..<nPoints
	let sincosarg = range.map {Double($0) * 1.0 * freq0 * dt}
	let (sine, cosine) = sincos(sincosarg) //Accelerate/Surge

	/* this gives normalization of the overlap sums */
  let sum_ee1 = zip(cosine, hanning).reduce(0, {$0 + $1.0*$1.0*$1.1})
	let sum_ee2 = zip(  sine, hanning).reduce(0, {$0 + $1.0*$1.0*$1.1})
	
  var sum_ef1 = zip(cosine, real).reduce(0, {$0 + $1.0*$1.1})
  sum_ef1 +=    zip(  sine, imag).reduce(0, {$0 + $1.0*$1.1})
  var sum_ef2 = zip(cosine, imag).reduce(0, {$0 + $1.0*$1.1})
  sum_ef2 -=    zip(  sine, real).reduce(0, {$0 + $1.0*$1.1})
  
	//let sum_ef1 = sum(mul(cosine, real)) + sum(mul(sine, imag))
	//let sum_ef2 = sum(mul(cosine, imag)) - sum(mul(sine, real))
	
  let sum1 = zip(real, imag).reduce(0.0, {$0 + $1.0*$1.0 + $1.1*$1.1})

  //(svesq(real)+svesq(imag)) //Accelerate/Surge
	//    NAFFData[i] -= (sum_ef1/sum_ee1*cosine[i] + sum_ef2/sum_ee2*sine[i])*hanning[i];
	
  let term1 = zip(cosine, sine).map {$0 * (sum_ef1/sum_ee1) + $1 * (sum_ef1/sum_ee1)}
  let term2 = zip(cosine, sine).map {$0 * (sum_ef2/sum_ee2) - $1 * (sum_ef2/sum_ee2)}

	//let term1 = add(vsmul(cosine, sum_ef1/sum_ee1),vsmul(sine, sum_ef1/sum_ee2)) //Accelerate/Surge
	//let term2 = sub(vsmul(cosine, sum_ef2/sum_ee2),vsmul(sine, sum_ef2/sum_ee2)) //Accelerate/Surge
  let term1win = zip(term1, hanning).map{$0 * $1}
  let term2win = zip(term2, hanning).map{$0 * $1}

  let newReal = zip(real, term1win).map { $0 - $1}
  let newImag = zip(imag, term2win).map { $0 - $1}
	
	let sum2 = zip(newReal, newImag).reduce(0.0, {$0 + $1.0*$1.0 + $1.1*$1.1}) //Accelerate/Surge
	
	let significance = (sum1>0.0) ? (sum2/sum1) : (-1.0)
	
	let freq0_2 = frequency / (Double.pi*2.0)
	let amplitude : Double = sqrt(pow(sum_ef1/sum_ee1,2.0)+pow(sum_ef2/sum_ee2,2.0))

	/* compute the phase and ensure that it is in the range [-PI, PI] */
	var phase = fmod(atan2(-sum_ef2/sum_ee2, sum_ef1/sum_ee1) +
		freq0_2*t0*Double.pi*2.0, Double.pi * 2.0);
	
	if (phase < -Double.pi){
		phase += Double.pi * 2.0
	}
	if (phase > Double.pi){
		phase -= Double.pi * 2.0
	}
	return (newReal, newImag,amplitude,phase,significance)
	
}


func CalculatePhaseAndAmplitudeFromFreqCplx(hanning: [Double], real: [Double], imag: [Double], nPoints: Int, dt: Double, frequency: Double, t0: Double)->(newReal:[Double], newImag: [Double], amplitude:Double, phase:Double,significance:Double){
    
    //var results = [Double](repeating: 0.0, count: data.count)
    
    assert (nPoints>1,"number of Points in CalculatePhaseAndAmplitudeFromFreq must be >1")
    assert (nPoints == hanning.count,"number of hanning window points must equal data points")
    
    let freq0 : Double = frequency  // * 2.0 * M_PI
    
    let range = 0..<nPoints
    let sincosarg = range.map {Double($0) * 1.0 * freq0 * dt}
    let (sine, cosine) = sincos(sincosarg) //Accelerate/Surge
    /* this gives normalization of the overlap sums */
    let sum_ee1 = zip(cosine, hanning).reduce(0.0, {$0 + $1.0*$1.0*$1.1})
    let sum_ee2 = zip(  sine, hanning).reduce(0.0, {$0 + $1.0*$1.0*$1.1})
    /* these are the overlap sums */
    let sum_ef1_r = zip(cosine, real).reduce(0.0, {$1.0*$1.1})
	  let sum_ef2_r = zip(  sine, real).reduce(0.0, {$1.0*$1.1})
	  let sum_ef1_i = zip(cosine, imag).reduce(0.0, {$1.0*$1.1})
	  let sum_ef2_i = zip(  sine, imag).reduce(0.0, {$1.0*$1.1})
	
    let sum1 = zip(real, imag).reduce(0.0, {$0 + $1.0*$1.0 + $1.1*$1.1})

  
  let term1 = zip(cosine, sine).map {$0 * (sum_ef1_r/sum_ee1) + $1 * (sum_ef2_r/sum_ee2)}
  let term2 = zip(cosine, sine).map {$0 * (sum_ef1_i/sum_ee1) + $1 * (sum_ef2_i/sum_ee2)}

  let term1win = zip(term1,hanning).map{$0 * $1}
  let term2win = zip(term2,hanning).map{$0 * $1}
  let newReal = zip(real, term1win).map{$0 - $1}
  let newImag = zip(imag, term2win).map{$0 - $1}
    
    let sum2 = zip(newReal, newImag).reduce(0.0, {$0 + $1.0*$1.0 + $1.1*$1.1})
    
    let significance = (sum1>0.0) ? (sum2/sum1) : (-1.0)
    
    let freq0_2 = frequency / (Double.pi*2.0)
    
    let amplitude : Double = sqrt((pow(sum_ef1_r/sum_ee1,2.0)+pow(sum_ef2_r/sum_ee2,2.0)+pow(sum_ef1_i/sum_ee1,2.0)+pow(sum_ef2_i/sum_ee2,2.0))/2.0)
    /* compute the phase and ensure that it is in the range [-PI, PI] */
    var phase = fmod(atan2(-sum_ef2_r/sum_ee2,sum_ef1_r/sum_ee1) + freq0_2*t0*Double.pi*2.0, Double.pi * 2.0)
    //var phase = fmod(atan2(-sum_ef2/sum_ee2, sum_ef1/sum_ee1) +
    //    freq0_2*t0*Double.pi*2.0, Double.pi * 2.0);
    //print(phase)
    if (phase < -Double.pi){
        phase += Double.pi * 2.0
    }
    if (phase > Double.pi){
        phase -= Double.pi * 2.0
    }
    return (newReal, newImag,amplitude,phase,significance)
    
}


/*
* optimize (e.g., minimize) a function of 1 variable using parabolic interpolation
*
* M. Borland, ANL/APS, 2001
*/
public enum parabolicOptResult{
	case Success(xBest:Double, yReturn:Double)
	case Failure(code:Int)
}


public func OneDParabolicOptimization(xGuess: Double, dx: Double, xLower: Double, xUpper: Double, function: (_ x:Double)->Double, maxCycles: Int, dxLimit: Double, tolerance: Double, maximize: Bool)->parabolicOptResult{
	let epsilon = Double.ulpOfOne
	
	
	var yReturn: Double
	let xGuess: Double = xGuess
	var dx = dx
	
	var maxFactor: Double = 1.0
	switch maximize{
	case true:
		maxFactor = -1.0
	case false:
		maxFactor = 1.0
	}
	
	var x0,x1,x2,xBest : Double
	var f0,f1,f2,fBest : Double
	x1 = 0.0
	f1 = 0.0
	
	x0 = xGuess
	f0 = maxFactor*function(x0)
	
	xBest = x0
	fBest = f0
	
	yReturn = maxFactor * f0
	
	/* find direction in which function decreases */
	var endCycle : Int = 0
	for cycle in 0..<(2*maxCycles){

		// let endCycle = cycle
		x1 = x0 + dx
		if (abs(x1-x0)<epsilon){
			//print(x0,x1,dx,epsilon)
		    hpNAFFLog("INFO: Parabolic search failure -1:x0,x1,dx,epsilon= \(x0),\(x1),\(dx),\(epsilon)", messageKind: .info)

			return .Failure(code: -1)
		}
		if ((x1>xUpper) || (x1<xLower)){
			return .Failure(code: -2)
		}
		f1 = maxFactor * function(x1)
		#if DEBUG
			//print("cycle = \(cycle), maxCycles=\(2*maxCycles), f1 = \(f1), fBest = \(fBest)")
		#endif
		if (f1<fBest){
			#if DEBUG
				//print("f1<fBest")
			#endif
			fBest = f1
			xBest = x1
		}
		if (f1<f0){
			#if DEBUG
				//print("f1<f0, breaking")
			#endif
			if (cycle == 2*maxCycles){
				if (fabs(dx)<dxLimit){
					return .Failure(code:-1)
				}
				return .Failure(code:-3)
			}
			break
		}
		dx *= (cycle%2==0 ? -1.0 : -0.5)
		if fabs(dx)<epsilon {
			break
		}
	}
	if (abs(x1-x0)<epsilon){
		return .Failure(code: -1)
	}
	#if DEBUG
		//print("Function decreases with dx=\(dx), f0=\(f0), f1=\(f1), cycle=\(endCycle)")
	#endif
	x2 = x1
	f2 = maxFactor * function(x2)
	/* take steps until passing through minimum */
	while (true){
		x2 = x1 + dx
		if ((x2>xUpper)||(x2<xLower)){
			return .Failure(code:-4)
		}
		f2 = maxFactor * function(x2)
		if (f2<fBest){
			fBest = f2
			xBest = x2
		}
		#if DEBUG
			//print("fBest = \(fBest), f1 = \(f1), f2 = \(f2)")
		#endif
		if (f2>f1){
			break
		}
		if (fabs(x1-x2)<epsilon){
			break
		}
		x0 = x1
		f0 = f1
		x1 = x2
		f1 = f2
	}
	if (x0>x2){
		/* arrange in increasing order */
		swap(&x0, &x2)
		swap(&f0, &f2)
	}
	/* now f0 > f1 and f2 > f1 */
	endCycle = 0
	for cycle in 0..<maxCycles {
		endCycle = cycle
		var x3,f3 : Double
		#if DEBUG
			//print("Cycle \(cycle):  f(\(x0))=\(f0),  f(\(x1))=\(f1),  f(\(x2))=\(f2)")
		#endif
		if ( ((x2-x0)<epsilon) || ((x2-x0)<dxLimit) || ((max(f2,f0)-f1)<tolerance)){
			break
		}
		/* try parabolic interpolation */
		let numer : Double = pow(x1-x0,2.0)*(f1-f2) - pow(x1-x2,2)*(f1-f0)
		let denom : Double = (x1-x0)*(f1-f2) - (x1-x2)*(f1-f0)
		x3 = x1 - numer/denom/2.0
		var failed : Bool = true
		let scale = x2-x0
		#if DEBUG
			//let debugf3 = isinf(x3) ? DBL_MAX : (maxFactor*function(x:x3))
			//print("parabolic parameters: x3 = \(x3), f3 = \( x3.isInfinite ? DBL_MAX : (maxFactor*function(x3)) ), scale=\(scale), x0=\(x0), x2=\(x2)")
		#endif
		
		if (!x3.isInfinite && x0<x3 && x3<x2 &&
			fabs(x3-x0)>epsilon*scale && fabs(x3-x1)>epsilon*scale &&
			fabs(x3-x2)>epsilon*scale) {
				/* evaluate at parabolic interpolation point */
				failed = false
				f3 = maxFactor*function(x3)
				if (f3<fBest) {
					fBest = f3
					xBest = x3
				}
				if (f3<f1) {
					/* replace point 1 */
					f1 = f3
					x1 = x3
				} else if ((f2>f0) && (f3<f2)) {
					/* replace point 2 */
					f2 = f3
					x2 = x3
					if (x2<x1) {
						swap(&x1, &x2)
						swap(&f1, &f2)
					}
				} else if (f2<f0 && f3<f0) {
					/* replace point 0 */
					f0 = f3
					x0 = x3
					if (x0>x1) {
						swap(&x0, &x1)
						swap(&f0, &f1)
					}
				} else{
					failed = true
				}
		}
		#if DEBUG
			//if (!failed) {print("Parabolic interpolation succeeded")}
		#endif
		if (failed == true){
			var right : Int = 0
			
			for other in 0...1{
				/* try dividing one of the intervals */
				if (fabs(x0-x1)<fabs(x1-x2)) {
					if (other == 0) {
						x3 = (x1+x2)/2.0
						right = 1
					} else {
						x3 = (x0+x1)/2.0
						right = 0
					}
				} else {
					if (other == 0) {
						x3 = (x0+x1)/2.0
						right = 0
					} else {
						x3 = (x1+x2)/2.0
						right = 1
					}
				}
				f3 = maxFactor*function(x3)
				if (f3<fBest) {
					fBest = f3
					xBest = x3
				}
				#if DEBUG
					//print("f3 = \(f3) at x3=\(x3)")
				#endif
				if (f3<f1) {
					/* replace point 1 */
					#if DEBUG
						//print("Replacing point 1")
					#endif
					f1 = f3
					x1 = x3
					break
				}
				if ( (right>0) && (f3<f2)) {
					/* replace point 2 */
					#if DEBUG
						//print("Replacing point 2")
					#endif
					f2 = f3
					x2 = x3
					if (x2<x1) {
						swap(&x1, &x2)
						swap(&f1, &f2)
					}
					break
				} else if ((right==0) && (f3<f0)) {
					/* replace point 0 */
					#if DEBUG
						//print("Replacing point 0")
					#endif
					f0 = f3
					x0 = x3
					if (x0>x1) {
						swap(&x0, &x1);
						swap(&f0, &f1);
					}
					break
				}
			}
			#if DEBUG
				//print("Sectioning succeeded")
			#endif
		}
	}
	#if DEBUG
		//print("Returning: x=\(xBest), y=\(maxFactor*fBest)")
	#endif
	yReturn = maxFactor * fBest
	
	return .Success(xBest: xBest, yReturn: yReturn)
	
}


public func performNAFF(data:[Double], cplxData:[Double]?=nil, dt: Double, nfreqs: Int, t0: Double, maxFrequencies: Int, fracRMSChangeLimit: Double, freqCycleLimit: Int, fracFreqAccuracyLimit: Double, lowerFreqLimit: Double,  upperFreqLimit: Double, complex: Bool = false) -> (frequencies:[Double],amplitudes:[Double],phases:[Double],significances:[Double]) {

	//var cplx = false
    //if let cplxData = cplxData {
    //    cplx = true
   // }
    
    
	let origCount = data.count
	assert(origCount > 1, "need at least 2 points in performNAFF")
	
	var NAFFData : [Double]
	var NAFFCplxData = [Double]()

	
	//Windowing
	var hanningwindow = hanning(points: origCount) //Accelerate/Surge
	//Remove mean
  let meanrealdata = data.reduce(0.0, {$0 + $1})/Double(data.count)
  //Surge.mean(data) //Accelerate/Surge
	let meanrealdatavec = [Double](repeating: meanrealdata, count: origCount)
  let windowedRealData0 = zip(data, meanrealdatavec).lazy.map{$0 - $1}
  let windowedRealData = zip(windowedRealData0, hanningwindow).map{$0 * $1}

	var meancplxdata : Double
	var windowedCplxData = [Double]()
	if let cplxData = cplxData, complex == true {
		meancplxdata = cplxData.reduce(0.0, {$0 + $1})/Double(cplxData.count)
		let meancplxdatavec = [Double](repeating: meancplxdata, count: origCount)
		//windowedCplxData = submul(cplxData, meancplxdatavec, hanningwindow) //Accelerate/Surge
    let windowedCplxData0 = zip(cplxData, meancplxdatavec).lazy.map{$0 - $1}
    windowedCplxData = zip(windowedCplxData0, hanningwindow).map{$0 * $1}
	}
    
	//Padding
	let paddedCount = nextDFTn(origCount)


	if (paddedCount == origCount){ // put window minimum in center of time series
		NAFFData = windowedRealData
		if (complex == true) {
			NAFFCplxData = windowedCplxData
		}
	}
	else
	{
		let paddingLength = paddedCount - origCount
		let padding = Array<Double>(repeating: 0.0, count: paddingLength)
		let paddingcenteridx = Int(ceil(Double(paddingLength)/2.0))
		// zero phase padding  http://dsp.stackexchange.com/questions/18938/merits-of-zero-phase-zero-padding/
		NAFFData = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), windowedRealData, Array(padding[paddingcenteridx..<paddingLength])].joined()))
		NAFFCplxData = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), windowedCplxData, Array(padding[paddingcenteridx..<paddingLength])].joined()))
		hanningwindow = Array(ContiguousArray([Array(padding[0..<paddingcenteridx]), hanningwindow, Array(padding[paddingcenteridx..<paddingLength])].joined()))

	}

	let count = NAFFData.count
	//let lowerFreqLimit = lowerFreqLimit ?? 0.0
	//let upperFreqLimit = upperFreqLimit ?? 1.0/(2.0*dt) // should be Nyquist frequency

	

	//_ = zip(NAFFData,NAFFCplxData).map{print("\($0.0)\t\($0.1)")}


	
	var frequencies = [Double](repeating:-1.0, count: maxFrequencies)
	var amplitudes = [Double](repeating:-1.0, count: maxFrequencies)
	var phases = [Double](repeating:-1.0, count: maxFrequencies)
	var significances = [Double](repeating:-1.0, count: maxFrequencies)
	//var residual = [Double](repeating:0.0, count: count)
	
	var NAFFdt : Double = dt //* Double(origCount-1)/Double(count-1)
	
	let freqSpacing = 1.0/(Double(count)*dt)
	
	
	//#warning("adjust for complex?")
	var rmsOrig: Double
	if (complex == true) {
	//rmsOrig = sqrt((svesq(NAFFData)+svesq(NAFFCplxData))/Double(count)) //sum of squares //Accelerate/Surge
		//rmsOrig = sqrt((svesq(NAFFData)+svesq(NAFFCplxData))/Double(count)) //sum of squares //Accelerate/Surge
    rmsOrig = sqrt(zip(NAFFData, NAFFCplxData).reduce(0.0, {$0 + $1.0*$1.0 + $1.1 * $1.1})/Double(count))
	} else
	{
    rmsOrig = sqrt(NAFFData.reduce(0.0, {$0 + $1*$1})/Double(count))
	}
	var rmsLast : Double = rmsOrig
	
	
	let FFTFreqs = (complex == true) ? count: count/2+1
	
	var freqsFound : Int = 0

	while(freqsFound < maxFrequencies){
		var magnitude2 : [Double]// = []
		if (complex == true) {
      let cplx = zip(NAFFData, NAFFData).map{Complex(real: $0, imag: $1)}
      magnitude2 = Fourier(cplx).map{$0.real*$0.real + $0.imag * $0.imag}
		} else {
      magnitude2 = Fourier(NAFFData).map{$0.real*$0.real + $0.imag * $0.imag}
		}
		var maxMag2 : Double = 0.0
		var iBest : Int = 0
		var ifreq_best = 0.0

		for i in 0..<FFTFreqs {
			var freq : Double
			if ((complex==true) && (i>(count/2))) {
				
					freq = -1.0*Double(count-i)*freqSpacing
			} else {
				freq = Double(i)*freqSpacing
			}
			
			if (magnitude2[i] > maxMag2){
				
				if ((lowerFreqLimit < upperFreqLimit) && (fabs(freq) < lowerFreqLimit) || (fabs(freq) > upperFreqLimit)) {
					continue
				}
				ifreq_best = freq
				iBest = i
				maxMag2 = magnitude2[i]
			}
		}
		if (iBest == 0){
			break
		}
		
		
		
		/*print(iBest)
		for i in 0..<count {
			print("\(i)\t\(magnitude2[i])")
		}*/
		//let wStart  = Double(iBest)*freqSpacing*Double.pi*2.0
		let wStart = ifreq_best*Double.pi*2.0
		frequencies[freqsFound] = wStart
		
		

		
		let dx = Double.pi*2.0*freqSpacing
		let xLower = (complex == true) ? -Double.pi/dt : 0.0
		let xUpper = Double.pi/dt
		func function(omega: Double) -> Double {
			if (complex == true) {
				return NAFFFunc2(NAFFData: NAFFData, NAFFCplxData: NAFFCplxData, NAFFdt: NAFFdt, omega: omega)
			} else {
				return NAFFFunc(NAFFData: NAFFData, NAFFdt: NAFFdt, omega: omega)
			}
		}
		
		let maxCycles = freqCycleLimit
		let dxLimit = Double.pi/NAFFdt * fracFreqAccuracyLimit
		let tolerance = 0.0
        //let scale = function(omega: wStart)

		//var res2 : parabolicOptResult
		// var newFreq : Double
		//var newAmplitude : Double
		
		let res = OneDParabolicOptimization(xGuess: frequencies[freqsFound], dx: dx, xLower: xLower, xUpper: xUpper, function: function, maxCycles: maxCycles, dxLimit: dxLimit, tolerance: tolerance, maximize: true)
		switch res{
		case .Success(let newFreq1, let amplitude1):
			//print("success: newFreq=\(newFreq1), amplitude=\(amplitude1)")
			frequencies[freqsFound] = newFreq1
			amplitudes[freqsFound] = amplitude1
			let res2 = OneDParabolicOptimization(xGuess: frequencies[freqsFound], dx: dx, xLower: xLower, xUpper: xUpper, function: function, maxCycles: maxCycles, dxLimit: dxLimit, tolerance: tolerance, maximize: true)
			
			switch res2{
			case .Success(let newFreq2, let amplitude2):
				//print("success: newFreq2=\(newFreq2), amplitude2=\(amplitude2)")
				frequencies[freqsFound] = newFreq2
				amplitudes[freqsFound] = amplitude2
				
			case .Failure(let code):
			    hpNAFFLog("INFO: Parabolic search returned code \(code)", messageKind: .warn)
//				print("Failure: code=\(code)")
				amplitudes[freqsFound] = -1.0
				frequencies[freqsFound] = -1.0
				break
			}
		case .Failure(let code):
//			print("Failure: code=\(code)")
		    hpNAFFLog("INFO: Parabolic search returned code \(code)", messageKind: .warn)
     		amplitudes[freqsFound] = -1.0
			frequencies[freqsFound] = -1.0
			break
		}
		
		
        if (complex == true) {
            (NAFFData, NAFFCplxData, amplitudes[freqsFound], phases[freqsFound], significances[freqsFound]) = CalculatePhaseAndAmplitudeFromFreqCplx(hanning: hanningwindow, real: NAFFData, imag: NAFFCplxData, nPoints: count, dt: NAFFdt, frequency: frequencies[freqsFound], t0: t0)
        } else {
            (NAFFData, amplitudes[freqsFound], phases[freqsFound], significances[freqsFound]) = CalculatePhaseAndAmplitudeFromFreq(hanning: hanningwindow, data: NAFFData, nPoints: count, dt: NAFFdt, frequency: frequencies[freqsFound], t0: t0)
        }
		
		
		frequencies[freqsFound] /= (Double.pi * 2.0)
		
		#if DEBUG
		//print("f=\(frequencies[freqsFound])  a=\(amplitudes[freqsFound])  p=\(phases[freqsFound])  s=\(significances[freqsFound])\n")
		#endif
		freqsFound += 1
		var rmsNow : Double = 0.0
		if (fracRMSChangeLimit > 0) {
			/* determine if residual is too small to bother with */
			if (complex == true) {
        rmsNow = sqrt(zip(NAFFData, NAFFCplxData).reduce(0.0, {$0 + $1.0*$1.0 + $1.1 * $1.1})/Double(count))
				//rmsNow = sqrt((/*Surge.*/svesq(NAFFData) +/*Surge.*/svesq(NAFFCplxData))/Double(count)) //Accelerate/Surge

			} else {
				rmsNow = sqrt(NAFFData.reduce(0.0, {$0 + $1*$1})/Double(count)) //sqrt(/*Surge.*/svesq(NAFFData)/Double(count)) //Accelerate/Surge

			}
			
			if ((rmsLast-rmsNow)/rmsOrig < fracRMSChangeLimit){
				break
			}
		}
		rmsLast = rmsNow;
	}
	
	frequencies = Array(frequencies[0..<freqsFound])
	amplitudes = Array(amplitudes[0..<freqsFound])
	phases = Array(phases[0..<freqsFound])
    significances = Array(significances[0..<freqsFound])
	
	
    if (paddedCount != origCount){
        //Adjust phases for switch to padded time series
        let paddingLength = paddedCount - origCount
        let initialpaddinglength = Int(ceil(Double(paddingLength)/2.0))
        var newphases = [Double]()
        newphases.reserveCapacity(phases.count)
        for (freq, phase) in zip(frequencies, phases) {
            let extraphase = freq * Double(initialpaddinglength)*dt
            let newphase = (extraphase - floor(extraphase))*2.0*Double.pi + phase
            newphases.append( fmod(newphase, 2.0*Double.pi))
        }
        phases = newphases
    }

	
	//add DC component
	frequencies = Array([[0.0], frequencies].joined())
	amplitudes = Array([[meanrealdata], amplitudes].joined())
	phases = Array([[0.0], phases].joined())
	significances = Array([[0.0], significances].joined())
	
	let result : ([Double],[Double],[Double],[Double]) = (frequencies, amplitudes, phases, significances)
	return result
}

public func nextDFTn(_ n: Int)-> Int {
  
  // calculate minimum necessary base 2 exponent
  let minbase2 = Int(ceil(log2(Double(n))))
  
  var allowedDFTCases = [Int]()
  
  //if base2 < 4, only 2**base2 allowed (for zrop, 5)
  
  for i in 1..<5 {
    allowedDFTCases.append(1<<i)
  }
  
  //if base2 >= 4, valid cases also include multiples of 1, 3, 5, 15
  if minbase2 >= 5 {
    for i in 5...minbase2 {
      let pow2 = 1<<i
      allowedDFTCases.append(pow2)
      allowedDFTCases.append(3*pow2)
      allowedDFTCases.append(5*pow2)
      allowedDFTCases.append(15*pow2)
    }
  }
  return allowedDFTCases.sorted(by: {$0 < $1}).filter({$0 >= n}).first!
}
