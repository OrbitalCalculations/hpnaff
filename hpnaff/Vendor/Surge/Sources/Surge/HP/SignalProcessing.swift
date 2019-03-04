//
//  SignalProcessing.swift
//  Surge
//
//  Created by Heiko Pälike on 21/05/2015.
//  Copyright (c) 2015 Heiko Pälike. All rights reserved.
//

import Foundation
import Accelerate
//http://stackoverflow.com/questions/21893481/hilbert-transform-analytical-signal-using-apples-accelerate-framework


public func hilbert(input:[Double])->([Double],[Double]) {
	
	assert(input.count>0,"need at least one point to calculate Hilbert")
	
	var hasOddLength : Bool
	var inputsignal : [Double] = input
	if (input.count % 2) != 0{
		inputsignal.append(0.0)
		hasOddLength = true
	}
	else
	{
		inputsignal = input
		hasOddLength = false
	}
	
	let length = inputsignal.count
	let halflength = length/2
	
	var even = [Double](repeating: 0.0, count:halflength)
	var odd = [Double](repeating: 0.0, count:halflength)
	var complex = DSPDoubleSplitComplex(realp:&even, imagp: &odd)
	
	
	UnsafePointer(inputsignal).withMemoryRebound(to: DSPDoubleComplex.self, capacity: halflength) { inputsignal in
		vDSP_ctozD(inputsignal, vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(halflength))
	}
	

	let forward = vDSP_DFT_zrop_CreateSetupD(nil, vDSP_Length(length), .FORWARD)
	let inverse = vDSP_DFT_zop_CreateSetupD(forward, vDSP_Length(length), .INVERSE)
	assert((forward != nil) && (inverse != nil),"FFT Setup Failed")

	

	var real = [Double](repeating: 0.0, count:length)
	var imag = [Double](repeating: 0.0, count: length)
	
	vDSP_DFT_ExecuteD(forward!, &even, &odd, &real, &imag)
    //  At this point, we have the forward real-to-complex DFT, which agrees with
    //  MATLAB up to a factor of two.  Since we want to double all but DC and NY
    //  as part of the Hilbert transform anyway, I'm not going to bother to
    //  unscale the rest of the frequencies -- they're already the values that
    //  we really want.  So we just need to move NY into the "right place",
    //  and scale DC and NY by 0.5.  The reflection frequencies are already
    //  zeroed out because the real-to-complex DFT only writes to the first n/2
    //  elements of real and imag.
	real[0] *= 0.5;	real[length/2] = 0.5 * imag[0];	imag[0] = 0.0
	
//	println("Stage 2:");	for i in 0..<length{println("\(real[i])+\(imag[i])i")}
	
	var hilbertreal = [Double](repeating: 0.0, count: length)
	var hilbertimag = [Double](repeating: 0.0, count: length)
    //  Now we have the completed hilbert transform up to a scale factor of n.
    //  We can unscale using vDSP_vsmulD.
	vDSP_DFT_ExecuteD(inverse!, &real, &imag, &hilbertreal, &hilbertimag)
    //  Now we have the completed hilbert transform up to a scale factor of n.
    //  We can unscale using vDSP_vsmulD.

	var scale = 1.0/Double(length)
	
    var hilbertrealinp = hilbertreal
    var hilbertimaginp = hilbertimag
	vDSP_vsmulD(&hilbertrealinp, 1, &scale, &hilbertreal, 1, vDSP_Length(length))
	vDSP_vsmulD(&hilbertimaginp, 1, &scale, &hilbertimag, 1, vDSP_Length(length))
	
	//println("Stage 2:");	for i in 0..<length{println("\(hilbertreal[i])+\(hilbertimag[i])i")}
	
	
	if hasOddLength == true {
		hilbertreal.removeLast()
		hilbertimag.removeLast()
	}
	
	vDSP_DFT_DestroySetupD(forward!)
	vDSP_DFT_DestroySetupD(inverse!)
	
	return (hilbertreal,hilbertimag)
    /*
    //Note that if you already have your DFT setups built and your storage allocated, the computation is extremely straightforward; the “real work” is just:
    
    vDSP_DFT_ExecuteD(forward, even, odd, real, imag);
    real[0] *= 0.5; real[n/2] = 0.5*imag[0]; imag[0] = 0.0;
    vDSP_DFT_ExecuteD(inverse, real, imag, hilbertreal, hilbertimag);
    double scale = 1.0/n; vDSP_vsmulD(hilbert, 1, &scale, hilbert, 1, 2*n);
*/
}

/*
#include <Accelerate/Accelerate.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
	const vDSP_Length n = 8;
	vDSP_DFT_SetupD forward = vDSP_DFT_zrop_CreateSetupD(NULL, n, vDSP_DFT_FORWARD);
	vDSP_DFT_SetupD inverse = vDSP_DFT_zop_CreateSetupD(forward, n, vDSP_DFT_INVERSE);
	//  Look like a typo?  The real-to-complex DFT takes its input separated into
	//  the even- and odd-indexed elements.  Since the real signal is [ 1, 2, 3, ... ],
	//  signal[0] is 1, signal[2] is 3, and so on for the even indices.
	double even[n/2] = { 1, 3, 5, 7 };
	double odd[n/2] = { 2, 4, 6, 8 };
	double real[n] = { 0 };
	double imag[n] = { 0 };
	vDSP_DFT_ExecuteD(forward, even, odd, real, imag);
	//  At this point, we have the forward real-to-complex DFT, which agrees with
	//  MATLAB up to a factor of two.  Since we want to double all but DC and NY
	//  as part of the Hilbert transform anyway, I'm not going to bother to
	//  unscale the rest of the frequencies -- they're already the values that
	//  we really want.  So we just need to move NY into the "right place",
	//  and scale DC and NY by 0.5.  The reflection frequencies are already
	//  zeroed out because the real-to-complex DFT only writes to the first n/2
	//  elements of real and imag.
	real[0] *= 0.5; real[n/2] = 0.5*imag[0]; imag[0] = 0.0;
	printf("Stage 2:\n");
	for (int i=0; i<n; ++i) printf("%f%+fi\n", real[i], imag[i]);
	
	double hilbert[2*n];
	double *hilbertreal = &hilbert[0];
	double *hilbertimag = &hilbert[n];
	vDSP_DFT_ExecuteD(inverse, real, imag, hilbertreal, hilbertimag);
	//  Now we have the completed hilbert transform up to a scale factor of n.
	//  We can unscale using vDSP_vsmulD.
	double scale = 1.0/n; vDSP_vsmulD(hilbert, 1, &scale, hilbert, 1, 2*n);
	printf("Stage 3:\n");
	for (int i=0; i<n; ++i) printf("%f%+fi\n", hilbertreal[i], hilbertimag[i]);
	vDSP_DFT_DestroySetupD(inverse);
	vDSP_DFT_DestroySetupD(forward);
	return 0;
}
Note that if you already have your DFT setups built and your storage allocated, the computation is extremely straightforward; the “real work” is just:

vDSP_DFT_ExecuteD(forward, even, odd, real, imag);
real[0] *= 0.5; real[n/2] = 0.5*imag[0]; imag[0] = 0.0;
vDSP_DFT_ExecuteD(inverse, real, imag, hilbertreal, hilbertimag);
double scale = 1.0/n; vDSP_vsmulD(hilbert, 1, &scale, hilbert, 1, 2*n);
*/


public enum BiQuadFilterType{
	case LowPass
	case HighPass
	case BandPass
	case Notch
	case Peak
	case LowShelf
	case HighShelf
}

public struct BiQuadFilter{
	let filterType : BiQuadFilterType
	let q : Double
	let fc : Double
	let peakGainDB : Double
	
	
	var a0 = 0.0
	var a1 = 0.0
	var a2 = 0.0
	var b1 = 0.0
	var b2 = 0.0
	
	var z1 = 0.0
	var z2 = 0.0
	
	public init(filterType: BiQuadFilterType, q: Double, fc: Double, peakGainDB: Double){
		self.filterType = filterType
		self.q = q
		self.fc = fc
		self.peakGainDB = peakGainDB
		self.a0 = 0.0
		(self.a0, self.a1, self.a2, self.b1, self.b2) = calcBiQuad()
	}
	
	
	public func process(input: [Double])->[Double]{
		let coefficients = [self.a0, self.a1, self.a2, self.b1, self.b2]
		let setup = vDSP_biquad_CreateSetupD(coefficients, vDSP_Length(coefficients.count))!

		var result = [Double](repeating:0.0, count:input.count)
		var delay = [Double](repeating:0.0, count:input.count)

		vDSP_biquadD(setup, &delay, input, 1, &result, 1, vDSP_Length(input.count))
		return(result)
	}
	
	func calcBiQuad()->(Double,Double,Double,Double,Double){
		let V = pow(10.0, fabs(self.peakGainDB) / 20.0)
		let K = tan(Double.pi * self.fc)
		
		var a0,a1,a2,b1,b2 : Double
		
		switch filterType{
		case .LowPass:
			let norm = 1.0 / (1.0 + K / self.q + K * K)
			a0 = K * K * norm
			a1 = 2.0 * a0
			a2 = a0
			b1 = 2.0 * (K * K - 1.0) * norm
			b2 = (1.0 - K/self.q + K*K) * norm
			return (a0,a1,a2,b1,b2)
		case .HighPass:
			let norm = 1.0 / (1.0 + K / self.q + K * K)
			a0 = 1.0 * norm
			a1 = -2.0 * a0
			a2 = a0
			b1 = 2.0 * (K * K - 1.0) * norm
			b2 = (1.0 - K/self.q + K*K) * norm
			return (a0,a1,a2,b1,b2)
		case .BandPass:
			let norm = 1.0 / (1.0 * K / self.q + K * K)
			a0 = K / self.q * norm
			a1 = 0.0
			a2 = -a0
			b1 = 2.0 * (K * K - 1.0) * norm
			b2 = (1.0 - K/self.q + K * K) * norm
			return (a0,a1,a2,b1,b2)
		case .Notch:
			let norm = 1.0 / (1.0 * K / self.q + K * K)
			a0 = (1.0 + K * K) * norm
			a1 = 2.0 * (K * K - 1.0) * norm
			a2 = a0
			b1 = a1
			b2 = (1.0 - K/self.q + K * K) * norm
			return (a0,a1,a2,b1,b2)
		case .Peak:
			if (self.peakGainDB >= 0.0){ // boost
				let norm = 1.0 / (1.0 + 1.0 / self.q * K + K * K)
				a0 = (1.0 + V/self.q * K + K * K) * norm
				a1 = 2.0 * (K * K - 1.0) * norm
				a2 = (1.0 - V / self.q * K + K * K) * norm
				b1 = a1
				b2 = (1.0 - 1.0/self.q * K + K * K) * norm
				return (a0,a1,a2,b1,b2)
			}
			else{ //cut
				let norm = 1.0 / (1.0 + V / self.q * K + K * K)
				a0 = (1.0 + 1.0/self.q * K + K * K) * norm
				a1 = 2.0 * (K * K - 1.0) * norm
				a2 = (1.0 - 1.0 / self.q * K + K * K) * norm
				b1 = a1
				b2 = (1.0 - V/self.q * K + K * K) * norm
				return (a0,a1,a2,b1,b2)
			}
		case .LowShelf:
			if (self.peakGainDB >= 0.0){ // boost
				let norm: Double = 1.0 / (1.0 + sqrt(2.0)*K + K*K)
				a0 = (1.0 + sqrt(2.0 * V) * K + V * K * K) * norm
				a1 = 2.0 * (V * K * K - 1.0) * norm
                a2 = (V * K * K) * norm
                a2 = a2 + (1.0 - sqrt(2.0 * V) * K)  * norm
                //a2 = (1.0 - sqrt(2.0 * V) * K + V * K * K) * norm
				b1 = 2.0 * (K * K - 1.0) * norm
				b2 = (1.0 - sqrt(2.0) * K + K * K) * norm
				return (a0,a1,a2,b1,b2)
			}
			else{ //cut
				var norm = 1.0
                norm = norm / (1.0 + sqrt(2.0*V) * K + V * K * K)
				a0 = (1.0 + sqrt(2.0) * K + K * K) * norm
				a1 = 2.0 * (K * K - 1.0) * norm
				a2 = (1.0 - sqrt(2.0) * K + K * K) * norm
				b1 = 2.0 * (V * K * K - 1.0) * norm
                b2 = norm
				b2 = b2 * (1.0 - (sqrt(2.0 * V) * K) + (V * K * K))
				return (a0,a1,a2,b1,b2)
			}
		case .HighShelf:
			if (self.peakGainDB >= 0.0){ // boost
				var norm : Double = 1.0
				norm /= (1.0 + sqrt(2.0)*K + K*K)
				a0 = (V + sqrt(2.0 * V) * K + K * K) * norm
				a1 = 2.0 * (K * K - V) * norm
				a2 = (V - sqrt(2.0 * V) * K + K * K) * norm
				b1 = 2.0 * (K * K - 1.0) * norm
				b2 = (1.0 - sqrt(2.0) * K + K * K) * norm
				return (a0,a1,a2,b1,b2)
			}
			else{ //cut
				let norm = 1.0 / (V + sqrt(2.0 * V) * K + K * K)
				a0 = (1.0 + sqrt(2.0) * K + K * K) * norm
				a1 = 2.0 * (K * K - 1.0) * norm
				a2 = (1.0 - sqrt(2.0) * K + K * K) * norm
				b1 = 2.0 * (K * K - V) * norm
				b2 = (V - sqrt(2.0 * V) * K + K * K) * norm
				return (a0,a1,a2,b1,b2)
			}
		}
	}
}

/*		switch (this->type) {

		case bq_type_highshelf:
			if (peakGain >= 0) {    // boost
				norm = 1 / (1 + sqrt(2) * K + K * K);
				a0 = (V + sqrt(2*V) * K + K * K) * norm;
				a1 = 2 * (K * K - V) * norm;
				a2 = (V - sqrt(2*V) * K + K * K) * norm;
				b1 = 2 * (K * K - 1) * norm;
				b2 = (1 - sqrt(2) * K + K * K) * norm;
			}
			else {    // cut
				norm = 1 / (V + sqrt(2*V) * K + K * K);
				a0 = (1 + sqrt(2) * K + K * K) * norm;
				a1 = 2 * (K * K - 1) * norm;
				a2 = (1 - sqrt(2) * K + K * K) * norm;
				b1 = 2 * (K * K - V) * norm;
				b2 = (V - sqrt(2*V) * K + K * K) * norm;
			}
			break;
		}
		
		return;
	}*/


/*

//http://www.earlevel.com/main/2012/11/26/biquad-c-source-code/
//
//  Biquad.h
//
//  Created by Nigel Redmon on 11/24/12
//  EarLevel Engineering: earlevel.com
//  Copyright 2012 Nigel Redmon
//
//  For a complete explanation of the Biquad code:
//  http://www.earlevel.com/main/2012/11/26/biquad-c-source-code/
//
//  License:
//
//  This source code is provided as is, without warranty.
//  You may copy and distribute verbatim copies of this document.
//  You may modify and use this source code to create binary code
//  for your own purposes, free or commercial.
//

#ifndef Biquad_h
#define Biquad_h

enum {
	bq_type_lowpass = 0,
	bq_type_highpass,
	bq_type_bandpass,
	bq_type_notch,
	bq_type_peak,
	bq_type_lowshelf,
	bq_type_highshelf
};

class Biquad {
	public:
	Biquad();
	Biquad(int type, double Fc, double Q, double peakGainDB);
	~Biquad();
	void setType(int type);
	void setQ(double Q);
	void setFc(double Fc);
	void setPeakGain(double peakGainDB);
	void setBiquad(int type, double Fc, double Q, double peakGain);
	float process(float in);
	
	protected:
	void calcBiquad(void);
	
	int type;
	double a0, a1, a2, b1, b2;
	double Fc, Q, peakGain;
	double z1, z2;
};

inline float Biquad::process(float in) {
	double out = in * a0 + z1;
	z1 = in * a1 + z2 - b1 * out;
	z2 = in * a2 - b2 * out;
	return out;
}

#endif // Biquad_h
//
//  Biquad.cpp
//
//  Created by Nigel Redmon on 11/24/12
//  EarLevel Engineering: earlevel.com
//  Copyright 2012 Nigel Redmon
//
//  For a complete explanation of the Biquad code:
//  http://www.earlevel.com/main/2012/11/26/biquad-c-source-code/
//
//  License:
//
//  This source code is provided as is, without warranty.
//  You may copy and distribute verbatim copies of this document.
//  You may modify and use this source code to create binary code
//  for your own purposes, free or commercial.
//

#include <math.h>
#include "Biquad.h"

Biquad::Biquad() {
	type = bq_type_lowpass;
	a0 = 1.0;
	a1 = a2 = b1 = b2 = 0.0;
	Fc = 0.50;
	Q = 0.707;
	peakGain = 0.0;
	z1 = z2 = 0.0;
}

Biquad::Biquad(int type, double Fc, double Q, double peakGainDB) {
	setBiquad(type, Fc, Q, peakGainDB);
	z1 = z2 = 0.0;
}

Biquad::~Biquad() {
}

void Biquad::setType(int type) {
	this->type = type;
	calcBiquad();
}

void Biquad::setQ(double Q) {
	this->Q = Q;
	calcBiquad();
}

void Biquad::setFc(double Fc) {
	this->Fc = Fc;
	calcBiquad();
}

void Biquad::setPeakGain(double peakGainDB) {
	this->peakGain = peakGainDB;
	calcBiquad();
}

void Biquad::setBiquad(int type, double Fc, double Q, double peakGainDB) {
	this->type = type;
	this->Q = Q;
	this->Fc = Fc;
	setPeakGain(peakGainDB);
}

void Biquad::calcBiquad(void) {
	double norm;
	double V = pow(10, fabs(peakGain) / 20.0);
	double K = tan(M_PI * Fc);
	switch (this->type) {
	case bq_type_lowpass:
		norm = 1 / (1 + K / Q + K * K);
		a0 = K * K * norm;
		a1 = 2 * a0;
		a2 = a0;
		b1 = 2 * (K * K - 1) * norm;
		b2 = (1 - K / Q + K * K) * norm;
		break;
		
	case bq_type_highpass:
		norm = 1 / (1 + K / Q + K * K);
		a0 = 1 * norm;
		a1 = -2 * a0;
		a2 = a0;
		b1 = 2 * (K * K - 1) * norm;
		b2 = (1 - K / Q + K * K) * norm;
		break;
		
	case bq_type_bandpass:
		norm = 1 / (1 + K / Q + K * K);
		a0 = K / Q * norm;
		a1 = 0;
		a2 = -a0;
		b1 = 2 * (K * K - 1) * norm;
		b2 = (1 - K / Q + K * K) * norm;
		break;
		
	case bq_type_notch:
		norm = 1 / (1 + K / Q + K * K);
		a0 = (1 + K * K) * norm;
		a1 = 2 * (K * K - 1) * norm;
		a2 = a0;
		b1 = a1;
		b2 = (1 - K / Q + K * K) * norm;
		break;
		
	case bq_type_peak:
		if (peakGain >= 0) {    // boost
			norm = 1 / (1 + 1/Q * K + K * K);
			a0 = (1 + V/Q * K + K * K) * norm;
			a1 = 2 * (K * K - 1) * norm;
			a2 = (1 - V/Q * K + K * K) * norm;
			b1 = a1;
			b2 = (1 - 1/Q * K + K * K) * norm;
		}
		else {    // cut
			norm = 1 / (1 + V/Q * K + K * K);
			a0 = (1 + 1/Q * K + K * K) * norm;
			a1 = 2 * (K * K - 1) * norm;
			a2 = (1 - 1/Q * K + K * K) * norm;
			b1 = a1;
			b2 = (1 - V/Q * K + K * K) * norm;
		}
		break;
	case bq_type_lowshelf:
		if (peakGain >= 0) {    // boost
			norm = 1 / (1 + sqrt(2) * K + K * K);
			a0 = (1 + sqrt(2*V) * K + V * K * K) * norm;
			a1 = 2 * (V * K * K - 1) * norm;
			a2 = (1 - sqrt(2*V) * K + V * K * K) * norm;
			b1 = 2 * (K * K - 1) * norm;
			b2 = (1 - sqrt(2) * K + K * K) * norm;
		}
		else {    // cut
			norm = 1 / (1 + sqrt(2*V) * K + V * K * K);
			a0 = (1 + sqrt(2) * K + K * K) * norm;
			a1 = 2 * (K * K - 1) * norm;
			a2 = (1 - sqrt(2) * K + K * K) * norm;
			b1 = 2 * (V * K * K - 1) * norm;
			b2 = (1 - sqrt(2*V) * K + V * K * K) * norm;
		}
		break;
	case bq_type_highshelf:
		if (peakGain >= 0) {    // boost
			norm = 1 / (1 + sqrt(2) * K + K * K);
			a0 = (V + sqrt(2*V) * K + K * K) * norm;
			a1 = 2 * (K * K - V) * norm;
			a2 = (V - sqrt(2*V) * K + K * K) * norm;
			b1 = 2 * (K * K - 1) * norm;
			b2 = (1 - sqrt(2) * K + K * K) * norm;
		}
		else {    // cut
			norm = 1 / (V + sqrt(2*V) * K + K * K);
			a0 = (1 + sqrt(2) * K + K * K) * norm;
			a1 = 2 * (K * K - 1) * norm;
			a2 = (1 - sqrt(2) * K + K * K) * norm;
			b1 = 2 * (K * K - V) * norm;
			b2 = (V - sqrt(2*V) * K + K * K) * norm;
		}
		break;
	}
	
	return;
}
*/


public func gaussian_filter(input: [Double], frequency: Double, bandwidth: Double, symmetric: Bool)->[Double]{


	return []

}

/*
//case k_sym_gaussianFilter:
{	size_t nn = 1;
if (GetCheckBoxValue( kFiltering_Notch_ChkBox ))	nn = 2;
double_vector1 fr(nn);	fr(nn) = freq;
double_vector1 bw(nn);	bw(nn) = bwidth;
double_vector1 c(nn);	c(nn)  = 1.0;
if (GetCheckBoxValue( kFiltering_Notch_ChkBox ))	{	fr(1) = 0;	bw(1) = 0;	c(1)  = 1.0;	c(nn) = -1.0;	};
if (popup == k_sym_gaussianFilter)
func = new MyMath::gaussian_filter_func( c, fr, bw, true );
else
func = new MyMath::gaussian_sum_func( c, fr, bw );
pxf = new double_vector1(3);	double_vector1& xf  = *pxf;
xf(1) = 0;		xf(2) = freq;		xf(3) = fc;
}
break;


class simplefunc : public doubleFunction {
	public:
	static	char*	TypeName( FourCharCode c )
	{	switch (c)
	{	case StairInterpFunctionCode:			return (char *)"Stair-mid";			break;
	case StairStartInterpFunctionCode:		return (char *)"Stair-start";		break;
	case StairEndInterpFunctionCode:		return (char *)"Stair-end";			break;
	case LinearInterpFunctionCode:			return (char *)"PiecewiseLinear";	break;
	case SplineInterpFunctionCode:			return (char *)"CubicSpline";		break;
	case PolynomialFunctionCode:			return (char *)"Polynomial";		break;
		}
		return	(char *)"";
	}
	protected:
	virtual double	func( const double x ) const {	return	ValueAt( x );	};
	public:
	
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) = 0;
	
	virtual double ValueAt( double t ) const = 0;
	virtual double DerivValueAt( double t ) const = 0;
	virtual double IntegrateBetween( double t1, double t2 ) const = 0;
	
	virtual FourCharCode code()	const = 0;
	
	virtual double MeanBetween( double t1, double t2 ) const
	{	return	( t1 == t2 ?  ValueAt( t1 ) : IntegrateBetween( t1, t2 )/(t2 - t1) );	};
};

class gaussian_sum_func : public simplefunc	{			//	sum of Ai exp( -( (f-fi)/wi )2 )
	protected:
	double_vector1 	coef;			//	Ai
	double_vector1 	freq;			//	fi
	double_vector1 	weight;			//	wi
	public:
	gaussian_sum_func( double_vector1& c, double_vector1& f, double_vector1& w ) :	coef(c), freq(f), weight(w)		{};
	gaussian_sum_func( double_vector1& f, double_vector1& w ) 	:	coef(1.0,f.length()), /*coef(f.length(),1.0), */freq(f), weight(w)		{};			//	coefs are set to 1
	
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new gaussian_sum_func( coef, freq, weight );	};
	
	virtual double ValueAt( double t ) const
		{	size_t	n = coef.length() + 1;
			double s = 0;
			while (--n>0)	(weight[n]!=0 ? s += coef[n] * exp( -sqr( (t-freq[n])/weight[n] ) ) : s += coef[n] );
			return	s;
	};
	virtual double DerivValueAt( double t ) const
		{	size_t	n = coef.length() + 1;
			double s = 0;
			while (--n>0)	(weight[n]!=0 ? s += coef[n] * 2*(freq[n]-t)/sqr(weight[n]) * exp( -sqr( (t-freq[n])/weight[n] ) ) : s += 0 );
			return	s;
	};
	virtual double IntegrateBetween( double t1, double t2 ) const
		{	size_t	n = coef.length() + 1;
			double s = 0;
			while (--n>0)	(weight[n]!=0 ? s += coef[n] * weight[n] * SqrtPiSur2 * (erf( (t2-freq[n])/weight[n] ) - erf( (t1-freq[n])/weight[n] )) : s += coef[n] * (t2-t1) );
			return	s;
	};
	
	virtual FourCharCode code()	const	{	return	GaussianSumFunctionCode;		};
};

class gaussian_filter_func : public gaussian_sum_func	{			//	idem 'gaussian_sum_func', with symetry and possible normalizations
	protected:
	void	enforce_symetry()					//	add the symetric ones
	{	size_t	n0 = coef.length();
		size_t	n = 2 * n0;
		coef.resize_but_keep( n );	freq.resize_but_keep( n );	weight.resize_but_keep( n );
		for (size_t	i=1; i<=n0; i++)
		{	coef[n0+i] = coef[i];	freq[n0+i] = -freq[i];	weight[n0+i] = weight[i];	}
	};
	void	enforce_max_to_one()
		{	size_t	n = coef.length();
			double	maxX = -1e222;
			for (size_t	i=1; i<=n; i++)
			{	double	x = ValueAt(freq[i]);	if (x>maxX)	maxX = x;	}		//	FAUX: le max. n'est pas 'exactement' a freq(i) des qu'il y a plusieurs gaussiennes.
			for (size_t	i=1; i<=n; i++)
				coef[i] /= maxX;
	};
	void	enforce_all_to_one()				//	recomputes new coef. so that result at freq(j) is one for each j.
		{	size_t	n = coef.length();			//	MAIS: le max. n'est pas 'exactement' a freq(i) des qu'il y a plusieurs gaussiennes.
			double_matrix1	m(n,n);
			double_vector1	v(n);
			int_vector1		indx(n);
			for (size_t	i=1; i<=n; i++)
			{	for (size_t	j=1; j<=n; j++)
				m(i,j) = exp( -sqr( (freq[j]-freq[i])/weight[i] ) );
				v(i) = 1;
			}
			Boolean even;
			ludcmp( m, indx, even );
			lubksb( m, indx, coef );
	};
	public:
	gaussian_filter_func( double_vector1& c, double_vector1& f, double_vector1& w, Boolean sym = false ) :	gaussian_sum_func( c, f, w )
	{	if (sym) enforce_symetry();
		enforce_max_to_one();
	};
	gaussian_filter_func( double_vector1& f, double_vector1& w, Boolean sym = false ) :	gaussian_sum_func( f, w )
	{	if (sym) enforce_symetry();
		enforce_all_to_one();
	};
	virtual simplefunc*	copySelf(const double_vector1&, const double_vector1&) 	{	return new gaussian_filter_func( coef, freq, weight );	};
	
	virtual FourCharCode code()	const	{	return	GaussianFilterFunctionCode;		};
};

*/


//http://hamiltonkibbe.com/finite-impulse-response-filters-using-apples-accelerate-framework-part-iii/
/*
/* filter signal x with filter h and store the result in output.
* output must be at least as long as x
*/

void
fft_fir_filter( const float *x,
	unsigned x_length,
	const float *h,
	unsigned h_length,
	float *output)
{
 
	// The length of the result from linear convolution is one less than the
	// sum of the lengths of the two inputs.
	unsigned result_length = x_length + h_length - 1;
	unsigned overlap_length = result_length - x_length;
	
 
	// Create buffer to store overflow across calls
	static float overflow[h_length - 1] = {0.0};
 
	// You need to implement next_pow2, it should return the first power of 2
	// greater than the argument passed to it.
	unsigned fft_length = next_pow2(result_length);
 
 
	// Create buffers to hold the zero-padded signal, filter kernel, and result.
	float h_padded[fft_length];
	float x_padded[fft_length];
	float temp_result[fft_length];
 
	// fill padded buffers with zeros
	float zero = 0.0;
	vDSP_vfill(&zero, h_padded, 1, fft_length);
	vDSP_vfill(&zero, x_padded, 1, fft_length);
 
	// Copy inputs into padded buffers
	cblas_scopy(h_length, h, 1, h_padded, 1);
	cblas_scopy(x_length, x, 1, x_padded, 1);
 
 
	// The Accelerate FFT needs an initialized setup structure. This, like much
	// of the above setup routine should be done outside of this function. I am
	// putting it here for ease of demonstration. This only needs to happen once.
	FFTSetup setup = vDSP_create_fftsetup(log2f((float)fft_length), FFT_RADIX2);
	
	// Create complex buffers for holding the Fourier Transforms of h and x
	// DSPSplitComplex holds pointers to two arrays of values, real, and imaginary.
	// Each array should hold at least fft_length/2 samples.
	DSPSplitComplex     h_DFT;
	DSPSplitComplex     x_DFT;
 
	// Create and assign the backing storage structures for each of these buffers.
	// In your actual implementation, these should be allocated elsewhere and
	// passed to this function along with the FFTSetup from above.
	float h_real[fft_length/2];
	float h_imag[fft_length/2];
	h_DFT.realp = h_real;
	h_DFT.imagp = h_imag;
	
	float x_real[fft_length/2];
	float x_imag[fft_length/2];
	x_DFT.realp = h_real;
	x_DFT.imagp = h_imag;
 
	
	// Here we calculate the FFT of the filter kernel. This should be done once
	// when you initialize your filter. As I mentioned previously, much of
	// this setup routine should be done outside of this function and saved.
 
	// Convert the real-valued filter kernel to split-complex form
	// and store it in our DFT array.
	vDSP_ctoz((DSPComplex*)h_padded, 2, &h_DFT, 1, (fft_length/2));
	
	// Do in-place FFT of filter kernel
	vDSP_fft_zrip(setup, &h_DFT, 1, log2f(fft_length), FFT_FORWARD);
 
	
	// Calculate FFT of input signal...
 
	// Convert the real-valued signal to split-complex form
	// and store it in our DFT array.
	vDSP_ctoz((DSPComplex*)x_padded, 2, &x_DFT, 1, (fft_length/2));
	
	// Do in-place FFT of the input signal
	vDSP_fft_zrip(setup, &x_DFT, 1, log2f(fft_length), FFT_FORWARD);
 
	
	// This gets a bit strange. The vDSP FFT stores the real value at nyquist in the
	// first element in the imaginary array. The first imaginary element is always
	// zero, so no information is lost by doing this. The only issue is that we are
	// going to use a complex vector multiply function from vDSP and it doesn't
	// handle this format very well. We calculate this multiplication ourselves and
	// add it into our result later.
	
	// We'll need this later
	float nyquist_mulitplied = h_DFT.imagp[0] * x_DFT.imagp[0];
 
	// Set the values in the arrays to 0 so they don't affect the multiplication
	h_DFT.imagp[0] = 0.0;
	x_DFT.imagp[0] = 0.0;
 
	// Complex multiply x_DFT by h_DFT and store the result in x_DFT
	vDSP_zvmul(&x_DFT, 1 &h_DFT, 1, &x_DFT,1, fft_length/2, 1);
 
	// Now we put our hand-calculated nyquist value back
	x_DFT.imagp[0] = nyquist_multiplied;
 
 
	// Do the inverse FFT of our result
	vDSP_fft_zrip(setup, &x_DFT, 1, log2f(fft_length), FFT_INVERSE);
 
	// And convert split-complex format to real-valued
	vDSP_ztoc(&x_DFT, 1, (DSPComplex*)temp_result, 2, fft_length/2);
 
	// Now we have to scale our result by 1/(4*fft_length)
	// (This is Apple's convention) to get our correct result.
	float scale = (1.0 / (4.0 * fft_length) );
	vDSP_vsmul(temp_result, 1, &scale, temp_result, 1, fft_length/2);
 
	
	// The rest of our overlap handling is the same as in our previous
	// implementation
 
	// Add the overlap from the previous run
	// use vDSP_vadd instead of loop
	vDSP_vadd(temp_result, overflow, temp_result, overlap_length);
	
	// Copy overlap into overlap buffer
	cblas_scopy(overlap_length, temp_result + x_length, 1, overflow, 1);
 
	// write the final result to the output. use BLAS copy instead of loop
	cblas_scopy(x_length, temp_result, 1, output, 1);
}
*/

public enum bandpassWindow {
	case Boxcar
	case Tukey(p: Double)
	case Gaussian(alpha: Double)
}

public enum bandpassPreprocess{
	case None
	case Demean
	case RemoveLinearTrend
}


func bandpassTaper(freqcount: Int, windowType: bandpassWindow)->[Double]{
	var taper = [Double](repeating: 1.0, count:freqcount)
	
	switch windowType{
	case .Boxcar:
		break
	case .Tukey(let p):
		let pbpts = Int(Double(freqcount) * p / 2.0)
		assert(pbpts>=2, "need to increase padding factor")
		//println(pbpts)
		
		for j in (0..<pbpts){
			let y = 2.0 * Double.pi * Double(j+1) / (p * Double(freqcount) + 1.0)
			taper[j] = 0.5 * (1.0 - cos(y))
			taper[freqcount - j - 1] = taper[j]
		}
	case .Gaussian(let alpha):
		// This portion of code is used if we want to apply a non-rectangular window
		// For a Gaussian window, based on Harris (1978), p.70. alpha is 1/std. dev.
		// (a measure of the width of the Fourier Transform of the Dirichlet kernel)
		// note that this function does not achieve zero, but approaches it for high alpha
		// here are some minima that correspond to different alpha:
		//  2.5 = 0.04711472; 3 = 0.01228416 ; 3.5 = 0.002508343.
		for j in (0..<freqcount){
			let x = (Double(j)+1.0) - (Double(freqcount)/2.0 + 0.5)
			let y = pow( alpha * x / (Double(freqcount)/2.0), 2.0)
			taper[j] = exp(-0.5 * y)
		}
	}

	return taper
}

public func bandpass(real_in: [Double], dt: Double, f_low: Double, f_high: Double, window: bandpassWindow = .Tukey(p: 0.4), preProcess:bandpassPreprocess = .RemoveLinearTrend, addMean: Bool = true, verbose: Bool = true)->[Double]  {
	let original_length = real_in.count
	
	// calculate Nyquist
	// let nyquist = 1.0 / (2.0 * dt)
	//calculate Rayleigh frequency
	// let rayleigh = 1.0 / (Double(original_length) * dt)
	//center bandpass frequency
	// let f_center = (f_high - f_low) / 2.0
	
	
	
	let dftLength = nextDFTn(original_length)
	let length = dftLength
	
	#if DEBUG
		//print("entering realfft, original length of input = \(original_length), using padded length \(dftLength)")
	#endif
	
	let paddingLength = dftLength - original_length
	let padding = [Double](repeating: 0.0, count: paddingLength)
	
	// let centeridx = Int(ceil(Double(original_length)/2.0))
	let paddingcenteridx = Int(ceil(Double(paddingLength)/2.0))

	
	//let inputsignal = [Double](real_in[0..<centeridx] + padding + real_in[centeridx..<original_length])
	

	
	
	// calculate mean value
	let meanValue = mean(real_in)
	
	// preprocess values (remove mean or linear trend)
	var inputsignal : [Double] = []
	
	switch preProcess{
	case .Demean:
        inputsignal = vsadd(real_in, -meanValue)
	case .RemoveLinearTrend:
		let (detrendedSignal, coefficients) = detrend(real_in, degree: 1)
		if (verbose == true) {print("function bandpass: Linear trend removed. m= \(coefficients.first!), intercept= \(coefficients.last!)")}
		inputsignal = detrendedSignal
	default:
		inputsignal = real_in
	}

	if (paddingLength > 0){
		inputsignal = Array([Array(padding[0..<paddingcenteridx]), inputsignal, Array(padding[paddingcenteridx..<paddingLength])].joined())
	}
	
	//inputsignal.extend(padding)
	
	
	let halflength = length/2
	
	
	
	let forward = vDSP_DFT_zrop_CreateSetupD(nil, vDSP_Length(length), .FORWARD)
	let inverse = vDSP_DFT_zrop_CreateSetupD(forward, vDSP_Length(length), .INVERSE)
	assert((forward != nil) && (inverse != nil),"FFT Setup Failed")
	
	

	var even = [Double](repeating: 0.0, count:halflength)
	var odd = [Double](repeating: 0.0, count:halflength)
	var complex = DSPDoubleSplitComplex(realp:&even, imagp: &odd)
	UnsafePointer(inputsignal).withMemoryRebound(to: DSPDoubleComplex.self, capacity: halflength) { inputsignal in
		vDSP_ctozD(inputsignal, vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(halflength))
	}
	
	//vDSP_ctozD(UnsafePointer<DSPDoubleComplex>(inputsignal), vDSP_Stride(2), &complex, vDSP_Stride(1), vDSP_Length(halflength))
	
	let fftoutputLength = halflength + 1 // have to leave one extra space to unpack vDSP DFT output
	
	var real = [Double](repeating: 0.0, count: fftoutputLength)
	var imag = [Double](repeating: 0.0, count: fftoutputLength)
	
	
	vDSP_DFT_ExecuteD(forward!, &even, &odd, &real, &imag)
	//  At this point, we have the forward real-to-complex DFT, which agrees with
	//  MATLAB up to a factor of two.  Since we want to double all but DC and NY
	//  as part of the Hilbert transform anyway, I'm not going to bother to
	//  unscale the rest of the frequencies -- they're already the values that
	//  we really want.  So we just need to move NY into the "right place",
	//  and scale DC and NY by 0.5.  The reflection frequencies are already
	//  zeroed out because the real-to-complex DFT only writes to the first n/2
	//  elements of real and imag.
	//real[0] *= 0.5;	real[length/2] = 0.5 * imag[0];	imag[0] = 0.0
	
	real[halflength] = imag[0]
	imag[0] = 0.0
	
	
	//var fftMagnitude = [Double](repeating: 0.0, count:halflength)
	
	//var fftoutput = DSPDoubleSplitComplex(realp: &real, imagp: &imag)

	//vDSP_zvmagsD(&fftoutput, 1, &fftMagnitude, 1, vDSP_Length(halflength)) // compute squared magnitude
	//fftMagnitude = vsmul(fftMagnitude, 1.0/4.0)
	
	//### new resulting frequency grid
	let nf = dftLength
	let df = 1.0 / (Double(nf) * dt)
	//let freq = (0...halflength).map({Double($0)*df})
	
	let mapped_f_low_idx = Int(ceil(f_low/df))
	let mapped_f_high_idx = Int(floor(f_high/df))
	let freqcount = (mapped_f_high_idx - mapped_f_low_idx) + 1

	var pbpts : Int
	switch window{
	case .Tukey(let p):
		pbpts = Int(Double(freqcount) * p / 2.0)
		assert(pbpts>=2, "need to increase padding factor")
	default:
		break
	}

	let taper = bandpassTaper(freqcount: freqcount, windowType: window)
	
	// apply window
	
	var filtered_real = [Double](repeating : 0.0, count:fftoutputLength)
	var filtered_imag = [Double](repeating : 0.0, count:fftoutputLength)
	
	filtered_real[mapped_f_low_idx...mapped_f_high_idx] = ArraySlice<Double>(mul([Double](real[mapped_f_low_idx...mapped_f_high_idx]),taper))
	filtered_imag[mapped_f_low_idx...mapped_f_high_idx] = ArraySlice<Double>(mul([Double](imag[mapped_f_low_idx...mapped_f_high_idx]),taper))
	
	// repack array
	
	filtered_imag[0] = real[halflength]
	filtered_real.removeLast()
	filtered_imag.removeLast()
	
	// reverse DFT
	var out_real = [Double](repeating: 0.0, count: halflength)
	var out_imag = [Double](repeating: 0.0, count: halflength)
	
	vDSP_DFT_ExecuteD(inverse!, &filtered_real, &filtered_imag, &out_real, &out_imag)
	
	var scale = 0.5/Double(length) //0.5 to account for vDSP zrop   impl = 2 * math
	
    var out_realinp = out_real
    var out_imaginp = out_imag
	vDSP_vsmulD(&out_realinp, 1, &scale, &out_real, 1, vDSP_Length(halflength))
	vDSP_vsmulD(&out_imaginp, 1, &scale, &out_imag, 1, vDSP_Length(halflength))
	
	var split_complex_out = DSPDoubleSplitComplex(realp:&out_real, imagp: &out_imag)
	
	var result = [Double](repeating: 0.0, count: length)
	
	UnsafeMutablePointer(mutating: result).withMemoryRebound(to: DSPDoubleComplex.self, capacity: halflength) { result in
		vDSP_ztocD(&split_complex_out, 1, result, 2, vDSP_Length(length))
	}

	//paddingcenteridx
	//paddingLength
	//original_length
	//result = [Double](result[(paddingcenteridx+1)...(paddingLength/2+original_length)])
	//result = [Double](result[0..<original_length])
	//plot(result, "result")
	//result.count

	
	if (addMean == true){
		result = vsadd(result, meanValue)
	}
	
	result = [Double](result[0..<original_length])

	
	vDSP_DFT_DestroySetupD(forward!)
	vDSP_DFT_DestroySetupD(inverse!)
	
	return result
}

/*
public func bandpass(t t: [Double],y inputData: [Double], paddingFactor: Int = 2, f_low: Double, f_high: Double, window: bandpassWindow = .Tukey(p:0.4), preProcess:bandpassPreprocess = .RemoveLinearTrend, addMean: Bool = true, verbose: Bool = true){
//re-implementation of astrochron bandpass function, file FUNCTION-bandpass_v13.R, Copyright (C) 2015 Stephen R. Meyers
	
	if (verbose == true) {print("\n----- BANDPASS FILTERING STRATIGRAPHIC SERIES-----\n")}
	
	// check input values for sanity
    assert((f_high > f_low), "function bandpass ERROR: parameter f_high should be larger than f_low")
	assert(t.count > 1,"function bandpass ERROR: need at least 2 points")
	assert(t.count == inputData.count, "function bandpass ERROR: number of data points in t needs to match those in y")
	
	let length = t.count
	let dt = t.last! - t.first!
	
	assert(dt >= 0, "function bandpass ERROR: data need to be sorted so that t is increasing")
	
    let dtest = sub([Double](t[1..<length]),y: [Double](t[0..<(length-1)]))
	let eps = max(dtest) - min(dtest)
	
	let epsm: Double = 1e-10 // let epsm = DBL_EPSILON
	assert(epsm > eps,"function bandpass ERROR: sampling interval is not uniform.")
	
	if (verbose == true) {print("function bandpass: Number of data points n= \(length), Sample Interval dt= \(dt)")}

	// calculate mean value
	let meanValue = mean(inputData)
	
	// preprocess values (remove mean or linear trend)
	var y: [Double] = []
    
	switch preProcess{
		case .None:
            y = inputData
		case .Demean:
            y = vsadd(inputData, c: -meanValue)
		case .RemoveLinearTrend:
			/*let (y, coefficients) = detrend(t, inputData, degree: 1)
			if (verbose == true) {println("function bandpass: Linear trend removed. m= \(coefficients.first!), intercept= \(coefficients.last!)")}*/
        break //HP_DEBUG

	}
	
	// calculate Nyquist
	let nyq = 1.0 / (2.0 * dt)
	//calculate Rayleigh frequency
	let ray = 1.0 / (Double(length) * dt)
	
	// pad with zeros if requested
    let newlength: Int
    if (((length * paddingFactor) % 2) != 0){ //check if odd: add another zero if we don't have an even number of data points, so Nyquist exists
         newlength = length * paddingFactor + 1
    }
    else
    {
        newlength = length * paddingFactor

    }
    let pad = y + [Double](count: (newlength-length), repeatedValue: 0.0)
    
    let nf = pad.count
    let df: Double = 1.0 / (Double(nf) * dt)
    var freq = [Double](count:nf, repeatedValue: 0.0)
    
    var f_highnew = f_high
    if (f_highnew > (nyq - df)){
        f_highnew = nyq - df
        print("function bandpass WARNING: resetting f_high to (nyq-df)=\(f_highnew) as it is too large")
    }
    
    var f_lownew = f_low
    if (f_lownew < df){
        f_lownew = df
        print("function bandpass WARNING: resetting f_low to df=\(f_lownew) as it is too small")
    }
    
    // center frequency of filter
    let fcent = (f_lownew + f_highnew)/2.0
    
    // take forward FFT
    let ft = fft(pad)
	
	
    
	/*
	### center of filter
	fcent=(flow+fhigh)/2
	
	### take fft
	ft <- fft(pad)
	# Locate real components for zero, nyquist, first neg freq., (zero-df)
	izero = 1
	nyqfreq = 0.5*nf + 1
	negfreq = 0.5*nf + 2
	minusdf = nf
	
	### set frequency index vector
	i <- seq(1,nf,by=1)
	### assign positive frequencies out to Nyquist
	freq <- df*(i[1:nyqfreq]-1)
	## assign negative frequencies
	f2 <- ( (nf/2) - (1:(minusdf-negfreq+1) ) ) * df * -1
	freq <- append(freq,f2)
	### caculate amplitude
	amp <- Mod(ft[1:nyqfreq])
	### put results into a data frame
	fft.out <- data.frame(cbind(freq[1:nyqfreq],amp))
	
	*/
	
	
	
}



/*
### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers


bandpass <- function (dat,padfac=2,flow=NULL,fhigh=NULL,win=0,alpha=3,p=0.25,demean=T,detrend=F,addmean=T,output=1,xmin=0,xmax=Nyq,genplot=T,verbose=T)
{
	
	if(verbose)   { cat("\n----- BANDPASS FILTERING STRATIGRAPHIC SERIES-----\n") }
	dat <- data.frame(dat)
	npts <- length(dat[,1])
	dt = dat[2,1]-dat[1,1]
	
	# error checking
	if(dt<0)
	{
		if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
		dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
		dt <- dat[2,1]-dat[1,1]
		npts <- length(dat[,1])
	}
	dtest <- dat[2:npts,1]-dat[1:(npts-1),1]
	epsm=1e-10
	if( (max(dtest)-min(dtest)) > epsm )
	{
		cat("\n**** ERROR: sampling interval is not uniform.\n")
		stop("**** TERMINATING NOW!")
	}
	
	d <- dat
	dt = d[2,1]-d[1,1]
	if(verbose)   { cat(" * Number of data points=", npts,"\n") }
	if(verbose)   { cat(" * Sample interval=", dt,"\n")}
	
	### remove mean
	dave <- colMeans(d[2])
	if (demean)
	{
		d[2] <- d[2] - dave
		if(verbose) { cat(" * Mean value removed=",dave,"\n") }
	}
	### use least-squares fit to remove linear trend
	if (detrend)
	{
		lm.0 <- lm(d[,2] ~ d[,1])
		d[2] <- d[2] - (lm.0$coeff[2]*d[1] + lm.0$coeff[1])
		if(verbose)  { cat(" * Linear trend removed. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n") }
	}
	
	### Calculate Nyquist
	Nyq <- 1/(2*dt)
	#### Calculate Rayleigh Frequency
	Ray <- 1/(npts*dt)
	### pad with zeros if desired (power of 2 not required!)
	### also convert from data frame to numeric
	pad <- as.numeric(d[,2])
	if(padfac>1) pad <- append( pad, rep(0,(npts*padfac-npts)) )
	
	### add another zero if we don't have an even number of data points, so Nyquist exists.
	if((npts*padfac)%%2 != 0) pad <- append(pad,0)
	
	### new resulting frequency grid
	nf = length(pad)
	df <- 1/(nf*dt)
	freq <- double(nf)
	
	if(is.null(flow)) flow <- df
	if(is.null(fhigh)) fhigh <- Nyq-df
	if (flow > fhigh)  {fh <- fhigh; fhigh <- flow; flow <- fh}
	if (fhigh> (Nyq-df))
	{
		fhigh <- (Nyq-df)
		if(verbose) cat(" * ERROR: maximum frequency too high, reset to Nyquist-df\n")
	}
	if (flow < df)
	{
		flow <- df
		if(verbose) cat(" * ERROR: minimum frequency too low, reset to df\n")
	}
	
	### center of filter
	fcent=(flow+fhigh)/2
	
	### take fft
	ft <- fft(pad)
	# Locate real components for zero, nyquist, first neg freq., (zero-df)
	izero = 1
	nyqfreq = 0.5*nf + 1
	negfreq = 0.5*nf + 2
	minusdf = nf
	
	### set frequency index vector
	i <- seq(1,nf,by=1)
	### assign positive frequencies out to Nyquist
	freq <- df*(i[1:nyqfreq]-1)
	## assign negative frequencies
	f2 <- ( (nf/2) - (1:(minusdf-negfreq+1) ) ) * df * -1
	freq <- append(freq,f2)
	### caculate amplitude
	amp <- Mod(ft[1:nyqfreq])
	### put results into a data frame
	fft.out <- data.frame(cbind(freq[1:nyqfreq],amp))
	
	### plot some of the results
	if(genplot)
	{
		par(mfrow=c(2,2))
		plot(d,type="l", col="blue",ylab="Value", xlab="Location", main="Stratigraphic Series",bty="n")
		plot(fft.out[,1],fft.out[,2], type="l",col="red",xlim=c(xmin,xmax),ylim=c(0,max(amp)),ylab="Amplitude", xlab="Frequency", main="Amplitude",bty="n")
		abline(v=flow,col="purple",lty=3)
		abline(v=fhigh,col="purple",lty=3)
	}
	
	### scaling factor for plotting of taper
	epsm=10^-13
	scaleT=max(amp)-min(amp)
	if(scaleT< epsm) scaleT = max(amp)
	
	### now bandpass
	bpts = 0
	# Apply rectangular window to zero and positive frequencies (first portion of the output)
	for (j in izero:nyqfreq) { if(freq[j] < flow || freq[j] > fhigh) {ft[j] = 0} else bpts=bpts+1 }
	# Apply rectangular window to negative frequencies (second portion of the output)
	for (j in negfreq:minusdf) { if(freq[j] < -1*fhigh || freq[j] > -1*flow) {ft[j] = 0} }
	if(verbose)
	{
		cat(" * Center of bandpass filter =",fcent,"\n")
		cat(" *",bpts,"pos/neg frequency pairs will be bandpassed\n")
	}
	
	### This portion of code is used if we want to apply a non-rectangular window
	### For a Gaussian window, based on Harris (1978), p.70. alpha is 1/std. dev.
	### (a measure of the width of the Fourier Transform of the Dirichlet kernel)
	### note that this function does not achieve zero, but approaches it for high alpha
	### here are some minima that correspond to different alpha:
	###  2.5 = 0.04711472; 3 = 0.01228416 ; 3.5 = 0.002508343.
	if(win == 1)
	{
		# make taper
		taper <- double(bpts)
		for (j in 1:bpts)
		{
			# center time series around x=0
			x= as.double(j) - ( ( as.double(bpts)/2 ) + 0.5)
			y= ( alpha * x / ( as.double(bpts)/2 ) )^2
			taper[j] = exp(-0.5*y)
		}
		# apply taper
		k = 1
		ii = 0
		for (j in izero:nyqfreq)
		{
			if(freq[j] >= flow && freq[j] <= fhigh)
			{
				ft[j] = ft[j]*taper[k]
				k=k+1
				if (ii == 0)
				{
					fhold <- freq[j]
					ii = ii +1
				}
			}
		}
		
		# now apply Gaussian taper to negative frequencies (note, this taper is symmetric, so the direction
		# doesn't matter)
		k = 1
		for (j in negfreq:minusdf) { if(freq[j] >= -1*fhigh && freq[j] <= -1*flow) {ft[j] = ft[j]*taper[k]; k=k+1} }
		# add window to spectrum plot
		i <- 1:bpts
		bpfreq <- fhold+df*(i-1)
		if(genplot)
		{
			lines(bpfreq,taper*scaleT,col="blue")
			points(fcent,max(amp),col="blue")
			mtext(round(fcent,digits=5),side=3,line=0,at=fcent,cex=0.5,font=4,col="blue")
		}
	}
	
	### For a simple cosine-tapered window. Based on Percival and Walden (1993), p. 209
	Pbpts=(bpts*p*.5)
	if(win == 2 && Pbpts >= 2)
	{
		# make taper
		taper <- double(bpts)
		taper[1:bpts] <- 1
		for (j in 1:Pbpts)
		{
			y= (2*pi*as.double(j) ) / ( (p*as.double(bpts) ) +1 )
			# First part of taper
			taper[j] = 0.5*(1-cos(y))
			# Last part of taper
			taper[bpts+1-j] = taper[j]
		}
		# apply taper
		k = 1
		ii = 0
		for (j in izero:nyqfreq)
		{
			if(freq[j] >= flow && freq[j] <= fhigh)
			{
				ft[j] = ft[j]*taper[k]
				k=k+1
				if (ii == 0)
				{
					fhold <- freq[j]
					ii = ii +1
				}
			}
		}
		# now apply cosine-taper to negative frequencies (note, this taper is symmetric, so the direction
		# doesn't matter)
		k = 1
		for (j in negfreq:minusdf) { if(freq[j] >= -1*fhigh && freq[j] <= -1*flow) {ft[j] = ft[j]*taper[k]; k=k+1} }
		# add window to spectrum plot
		i <- 1:bpts
		bpfreq <- fhold+df*(i-1)
		if(genplot)
		{
			lines(bpfreq,taper*scaleT,col="blue")
			points(fcent,max(amp),col="blue")
			mtext(round(fcent,digits=5),side=3,line=0,at=fcent,cex=0.5,font=4,col="blue")
		}
	}
	
	# if fewer than 4 tapered points (2 on each side of window), to not apply taper
	if(win ==2 && Pbpts < 2)
	{
		cat("\n**** WARNING: Too few data points to apply cosine-tapered window.\n")
		cat("              Will use rectangular window. Try increasing padfac,\n")
		cat("              and/or increasing p.\n")
	}
	
	# inverse FFT
	###  convert to real number and normalize iFFT output by length of record, suppress warning
	###   about discarding imaginary component (it is zero!).
	ifft <- suppressWarnings(as.double(fft(ft,inverse=TRUE)/nf))
	### isolate prepadded portion
	bp <- ifft[1:npts]
	d[2] <- bp
	d3<- data.frame( cbind(d[1],(bp + dave)) )
	colnames(d)[2] <- colnames(dat[2])
	colnames(d3)[2] <- colnames(dat[2])
	
	if(genplot)
	{
		plot(d, type="l",col="red", ylab="Value", xlab="Location", main="Bandpassed Signal",bty="n")
		plot(dat, type="l",col="blue", ylim=c(min(dat[,2],d3[,2]),max(dat[,2],d3[,2])), ylab="Value", xlab="Location", main="Comparison",bty="n")
		lines(d3, col="red")
	}
	
	if(output == 1)
	{
		if(!addmean) {return(d)}
		if(addmean) {return(d3)}
	}
	
	if(output == 2) return(data.frame(cbind(bpfreq,taper)))
	
	#### END function bandpass
}
*/



/*
For a complex DFT, the way the vDSP_DFT_Execute function works is fairly straightforward: it works on real/imaginary input/outputs specified by __vDSP_Ir, __vDSP_Ii, __vDSP_Or, and __vDSP_Oi, in the way indicated by their names and the transform is either forward or inverse according to the flag specified in the setup. For a real-to-complex DFT, there are some complications: forward means real-to-complex, inverse is complex-to-real. In either direction, the data for the real side of the transform is stored alternately in the real and imaginary arrays (even-indexed elements in the real array, odd-indexed in the imaginary array).

If the setup was created with vDSP_DFT_zop_CreateSetup and N was passed for the __vDSP_Length parameter, each vector should contain length elements. If the setup was created with vDSP_DFT_zrop_CreateSetup and N was passed for the __vDSP_Length, each vector should contain N/2 elements.

When the vDSP_DFT_Execute function is called with a setup returned from the vDSP_DFT_zop_CreateSetup function, it calculates:

for 0 <= k < N,
H[k] = sum(1**(S * j*k/N) * h[j], 0 <= j < N),
where:

N is the length given in the setup.

h is the array of complex numbers specified by Ir and Ii when the vDSP_DFT_Execute function is called:

for 0 <= j < N,
h[j] = Ir[j] + i * Ii[j];
H is the array of complex numbers specified by Or and Oi when the vDSP_DFT_Execute function returns:

for 0 <= k < N,
H[k] = Or[k] + i * Oi[k];
S is -1 if Direction is vDSP_DFT_FORWARD and +1 if Direction is vDSP_DFT_INVERSE; and

1**x is e**(2*pi*i*x).

When the vDSP_DFT_Execute function is called with a setup returned from the DSP_DFT_zrop_CreateSetup function, it calculates:

for 0 <= k < N,
H[k] = C * sum(1**(S * j*k/N) * h[j], 0 <= j < N),
where:

N is the length given in the setup.

h is the array of numbers specified by Ir and Ii when the vDSP_DFT_Execute is called (see "Data Layout" below).

H is the array of numbers specified by Or and Oi when the vDSP_DFT_Execute function returns (see "Data Layout", below.)

S is -1 if Direction is vDSP_DFT_FORWARD and +1 if Direction is vDSP_DFT_INVERSE; and

1**x is e**(2*pi*i*x).

Data Layout:

If Direction is vDSP_DFT_FORWARD, then:

h is an array of real numbers, with its even-index elements stored in Ir and its odd-index elements stored in Ii:

for 0 <= j < N/2,
h[2*j+0] = Ir[j], and
h[2*j+1] = Ii[j]
H is an array of complex numbers, stored in Or and Oi:

H[0  ] = Or[0].  // (H[0  ] is pure real.)
H[N/2] = Oi[0].  // (H[N/2] is pure real.)
for 1 < k < N/2
H[k] = Or[k] + i * Oi[k]
For N/2 < k < N, H[k] is not explicitly stored in memory but is known because it necessarily equals the conjugate of H[N-k], which is stored as described above.

If Direction is vDSP_DFT_Inverse, then the layouts of the input and output arrays are swapped. Ir and Ii describe an input array with complex elements laid out as for Or and Oi. When the vDSP_DFT_Execute function returns, Or and Oi contain a pure real array, with its even-index elements stored in Or and its odd-index elements in Oi.

*/

*/
