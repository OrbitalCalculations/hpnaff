// swift-tools-version:5.5
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "hpNaff",
    platforms: [.macOS(.v11)],
    dependencies: [
      //.package(url: "https://github.com/Jounce/Surge.git", .upToNextMajor(from: "2.3.2")),
      .package(url: "https://github.com/onevcat/Rainbow", .upToNextMajor(from: "4.0.0")),
      .package(url: "https://github.com/SusanDoggie/Doggie.git", .upToNextMajor(from: "6.6.1")),
      .package(url: "https://github.com/apple/swift-numerics", from: "1.0.0"),

        // Dependencies declare other packages that this package depends on.
        // .package(url: /* package url */, from: "1.0.0"),
    ],
    targets: [
        // Targets are the basic building blocks of a package. A target can define a module or a test suite.
        // Targets can depend on other targets in this package, and on products in packages this package depends on.
        .executableTarget(
            name: "hpNaff",
            dependencies: [
              "Rainbow", /*"Surge",*/
              .product(name: "DoggieMath", package: "Doggie"),
              .product(name: "Numerics", package: "swift-numerics"),

            ]),
        .testTarget(
            name: "hpNaffTests",
            dependencies: ["hpNaff"]),
    ]
)
