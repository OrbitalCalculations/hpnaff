//
//  Utilities.swift
//  hpnaff
//
//  Created by Heiko Pälike on 25/05/2018.
//  Copyright © 2018 Heiko Pälike. All rights reserved.
//

import Foundation

public enum logMessageType {
    /** About text */
    case standard

    /** An error message */
    case error

    /** A warn message */
    case warn

    /** An info message */
    case info
  }


public func hpNAFFLog(_ str: String, messageKind: logMessageType = .standard, verbosity: Int = 0 ) {
	switch messageKind {
	case .error:
	  NSLog(str.red)
	case .info:
	  if verbosity > 1 {
		NSLog(str.lightBlue)
	  }
    case .warn:
	  if verbosity > 0 {
		NSLog(str.lightRed)
	  }
	default:
	  NSLog(str.lightBlack)
	}
}

