// swift-tools-version:5.5

import PackageDescription

let package = Package(
    name: "hpNaff",
    platforms: [.macOS(.v11)],
    products: [
      .executable(name: "hpnaff", targets: ["hpNaff"])
    ],
    dependencies: [
      .package(url: "https://github.com/onevcat/Rainbow", .upToNextMajor(from: "4.0.0")),
      .package(url: "https://github.com/apple/swift-numerics", from: "1.0.0"),
      .package(url: "https://github.com/apple/swift-algorithms", from: "1.0.0"),
      .package(url: "https://github.com/apple/swift-collections", from: "1.0.0"),
    ],
    targets: [
        // Targets are the basic building blocks of a package. A target can define a module or a test suite.
        // Targets can depend on other targets in this package, and on products in packages this package depends on.
        .executableTarget(
            name: "hpNaff",
            dependencies: [
              "Rainbow",
              .product(name: "Numerics", package: "swift-numerics"),
              .product(name: "Algorithms", package: "swift-algorithms"),
              .product(name: "Collections", package: "swift-collections"),
            ],
            resources: [
              .process("Resources")
            ]
        ),
        .testTarget(
            name: "hpNaffTests",
            dependencies: ["hpNaff"]),
    ]
)
