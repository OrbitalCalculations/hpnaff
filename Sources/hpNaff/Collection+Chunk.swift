//
//  Collection+Chunk.swift
//  hpnaff
//
//  Created by Heiko Pälike on 25/05/2018.
//  Copyright © 2018 Heiko Pälike. All rights reserved.
//

import Foundation

public extension Collection {

  func chunked(by distance: Int) -> [[Element]] {
    var result: [[Element]] = []
    var batch: [Element] = []
    for element in self {
      batch.append(element)
      if batch.count == distance {
        result.append(batch)
        batch = []
      }
    }
    if !batch.isEmpty {
      result.append(batch)
    }
    return result
  }
	
	func chunked(by distance: Int, windowOffset: Int) -> [[Element]] {
        assert((windowOffset>0) && (windowOffset<=distance))
        var result: [[Element]] = []
        var batch: [Element] = []

  	var offset = 0
		var i = self.startIndex
		
		while self.index(i, offsetBy: offset) < self.endIndex {
			let index = self.index(i, offsetBy: offset)
			batch.append(self[index])
      i = self.index(after: i)
			
			if batch.count == distance {
        result.append(batch)
        batch = []
        offset += windowOffset
        i = self.startIndex
      }
		}
		
		if !batch.isEmpty {
      result.append(batch)
    }
    return result
  }
}
