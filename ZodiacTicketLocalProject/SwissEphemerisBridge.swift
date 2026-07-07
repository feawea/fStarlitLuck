import Foundation

enum SwissEphemerisError: Error {
    case ephePathMissing
    case calculationFailed(Int32)
    case timeConversionFailed
}

enum AstroBody: CaseIterable {
    case sun
    case moon
    case mercury
    case venus
    case mars
    case jupiter
    case saturn
    case uranus
    case neptune
    case pluto

    var code: Int32 {
        switch self {
        case .sun: return 0
        case .moon: return 1
        case .mercury: return 2
        case .venus: return 3
        case .mars: return 4
        case .jupiter: return 5
        case .saturn: return 6
        case .uranus: return 7
        case .neptune: return 8
        case .pluto: return 9
        }
    }
}

enum HouseSystem: String, CaseIterable, Codable, Identifiable {
    case placidus = "P"
    case koch = "K"
    case equal = "E"

    var id: String { rawValue }

    var swissCode: Int32 {
        Int32(rawValue.utf8.first ?? 80)
    }
}

struct PlanetPosition {
    let body: Int32
    let longitude: Double
    let latitude: Double
    let distance: Double
    let speedLongitude: Double

    var signIndex: Int {
        let normalized = AstroMath.normalizedDegrees(longitude)
        return Int(floor(normalized / 30.0)) % 12
    }

    var degreeInSign: Double {
        let normalized = AstroMath.normalizedDegrees(longitude)
        return normalized.truncatingRemainder(dividingBy: 30.0)
    }
}

struct HouseSet {
    let system: String
    let cusps: [Double]
    let ascendant: Double
    let mc: Double
}

struct HousePositions {
    let system: HouseSystem
    let cusps: [Double]
    let ascendant: Double
    let mc: Double
    let lstDegrees: Double
}

struct SwissNatalChart {
    let julianDayUT: Double
    let localSiderealTime: Double
    let planets: [PlanetPosition]
    let houses: HousePositions
}

struct BirthContext {
    let date: Date
    let latitude: Double
    let longitude: Double
    let timeZone: TimeZone
}

final class SwissEphemerisService {
    static let shared = SwissEphemerisService()
    private static let apiLock: NSRecursiveLock = {
        let lock = NSRecursiveLock()
        lock.name = "com.fdivination.swiss-ephemeris"
        return lock
    }()

    private var epheConfigured = false

    private init() {}

    private static func withExclusiveAccess<T>(_ body: () throws -> T) rethrows -> T {
        apiLock.lock()
        defer { apiLock.unlock() }
        return try body()
    }

    func configureEphemerisPath() throws {
        try Self.withExclusiveAccess {
            try configureEphemerisPathLocked()
        }
    }

    private func configureEphemerisPathLocked() throws {
        guard let path = Self.locateEphePath() else { throw SwissEphemerisError.ephePathMissing }
        swe_set_ephe_path(path)
        epheConfigured = true
    }

    private func ensureConfiguredLocked() throws {
        if !epheConfigured {
            try configureEphemerisPathLocked()
        }
    }

    func computeNatalChart(context: BirthContext, houseSystem: HouseSystem) throws -> SwissNatalChart {
        try Self.withExclusiveAccess {
            try ensureConfiguredLocked()
            let jdUT = try Self.julianDayUTLocked(for: context)
            let planets = try computePlanetsLocked(jdUT: jdUT, flags: SEFLG_SWIEPH | SEFLG_SPEED)
            let houseSet = try computeHousesLocked(
                jdUT: jdUT,
                latitude: context.latitude,
                longitude: context.longitude,
                system: houseSystem
            )
            let lst = Self.localSiderealTimeLocked(jdUT: jdUT, longitude: context.longitude)
            let houses = HousePositions(
                system: houseSystem,
                cusps: houseSet.cusps,
                ascendant: houseSet.ascendant,
                mc: houseSet.mc,
                lstDegrees: lst
            )
            return SwissNatalChart(julianDayUT: jdUT, localSiderealTime: lst, planets: planets, houses: houses)
        }
    }

    func computePlanets(jdUT: Double, flags: Int32 = SEFLG_SWIEPH | SEFLG_SPEED) throws -> [PlanetPosition] {
        try Self.withExclusiveAccess {
            try ensureConfiguredLocked()
            return try computePlanetsLocked(jdUT: jdUT, flags: flags)
        }
    }

    func computeHouses(jdUT: Double, latitude: Double, longitude: Double, system: HouseSystem) throws -> HouseSet {
        try Self.withExclusiveAccess {
            try ensureConfiguredLocked()
            return try computeHousesLocked(jdUT: jdUT, latitude: latitude, longitude: longitude, system: system)
        }
    }

    func computeAllHouses(jdUT: Double, latitude: Double, longitude: Double) throws -> [HouseSet] {
        try Self.withExclusiveAccess {
            try ensureConfiguredLocked()
            return [
                try computeHousesLocked(jdUT: jdUT, latitude: latitude, longitude: longitude, system: .placidus),
                try computeHousesLocked(jdUT: jdUT, latitude: latitude, longitude: longitude, system: .koch),
                try computeHousesLocked(jdUT: jdUT, latitude: latitude, longitude: longitude, system: .equal)
            ]
        }
    }

    func calcPlanet(jdUT: Double, planet: Planet) throws -> (lon: Double, speed: Double) {
        try Self.withExclusiveAccess {
            try ensureConfiguredLocked()
            var xx = [Double](repeating: 0, count: 6)
            var serr = [Int8](repeating: 0, count: 256)
            let flags = SEFLG_SWIEPH | SEFLG_SPEED
            let rc = swe_calc_ut(jdUT, planet.swissCode, flags, &xx, &serr)
            if rc < 0 { throw SwissEphemerisError.calculationFailed(rc) }
            return (lon: xx[0], speed: xx[3])
        }
    }

    func calcHouseDetails(
        jdUT: Double,
        latitude: Double,
        longitude: Double,
        system: HouseSystem
    ) throws -> (cusps: [Double], asc: Double, mc: Double, armc: Double, trueObliquity: Double) {
        try Self.withExclusiveAccess {
            try ensureConfiguredLocked()
            var cusps = [Double](repeating: 0, count: 13)
            var ascmc = [Double](repeating: 0, count: 10)
            var eclNut = [Double](repeating: 0, count: 6)
            var serr = [Int8](repeating: 0, count: 256)
            let rc = swe_houses_ex(jdUT, SEFLG_SWIEPH, latitude, longitude, system.swissCode, &cusps, &ascmc)
            if rc < 0 { throw SwissEphemerisError.calculationFailed(rc) }
            let eclRc = swe_calc_ut(jdUT, SE_ECL_NUT, SEFLG_SWIEPH, &eclNut, &serr)
            if eclRc < 0 { throw SwissEphemerisError.calculationFailed(eclRc) }
            return (
                cusps: cusps,
                asc: ascmc[Int(SE_ASC)],
                mc: ascmc[Int(SE_MC)],
                armc: ascmc[Int(SE_ARMC)],
                trueObliquity: eclNut[0]
            )
        }
    }

    static func julianDayUT(for context: BirthContext) throws -> Double {
        try withExclusiveAccess {
            try julianDayUTLocked(for: context)
        }
    }

    private static func julianDayUTLocked(for context: BirthContext) throws -> Double {
        var calendar = Calendar(identifier: .gregorian)
        calendar.timeZone = context.timeZone
        let comps = calendar.dateComponents([.year, .month, .day, .hour, .minute, .second], from: context.date)
        guard let year = comps.year, let month = comps.month, let day = comps.day,
              let hour = comps.hour, let minute = comps.minute else {
            throw SwissEphemerisError.timeConversionFailed
        }
        let seconds = Double(comps.second ?? 0)
        let offsetHours = Double(context.timeZone.secondsFromGMT(for: context.date)) / 3600.0

        var utcYear: Int32 = 0
        var utcMonth: Int32 = 0
        var utcDay: Int32 = 0
        var utcHour: Int32 = 0
        var utcMinute: Int32 = 0
        var utcSecond: Double = 0

        swe_utc_time_zone(
            Int32(year), Int32(month), Int32(day),
            Int32(hour), Int32(minute), seconds,
            offsetHours,
            &utcYear, &utcMonth, &utcDay,
            &utcHour, &utcMinute, &utcSecond
        )

        var jdResults = [Double](repeating: 0, count: 2)
        var serr = [Int8](repeating: 0, count: 256)
        let rc = swe_utc_to_jd(utcYear, utcMonth, utcDay, utcHour, utcMinute, utcSecond, Int32(SE_GREG_CAL), &jdResults, &serr)
        if rc < 0 { throw SwissEphemerisError.calculationFailed(rc) }

        return jdResults[1]
    }

    static func localSiderealTime(jdUT: Double, longitude: Double) -> Double {
        withExclusiveAccess {
            localSiderealTimeLocked(jdUT: jdUT, longitude: longitude)
        }
    }

    private static func localSiderealTimeLocked(jdUT: Double, longitude: Double) -> Double {

        let gstHours = swe_sidtime(jdUT)
        let lstHours = AstroMath.normalizedHours(gstHours + longitude / 15.0)
        return lstHours * 15.0
    }

    private func computePlanetsLocked(jdUT: Double, flags: Int32) throws -> [PlanetPosition] {
        var positions: [PlanetPosition] = []
        let bodies: [AstroBody] = AstroBody.allCases
        var serr = [Int8](repeating: 0, count: 256)

        for body in bodies {
            var xx = [Double](repeating: 0, count: 6)
            let rc = swe_calc_ut(jdUT, body.code, flags, &xx, &serr)
            if rc < 0 { throw SwissEphemerisError.calculationFailed(rc) }
            let pos = PlanetPosition(
                body: body.code,
                longitude: xx[0],
                latitude: xx[1],
                distance: xx[2],
                speedLongitude: xx[3]
            )
            positions.append(pos)
        }
        return positions
    }

    private func computeHousesLocked(
        jdUT: Double,
        latitude: Double,
        longitude: Double,
        system: HouseSystem
    ) throws -> HouseSet {
        var cusps = [Double](repeating: 0, count: 13)
        var ascmc = [Double](repeating: 0, count: 10)
        let rc = swe_houses_ex(jdUT, SEFLG_SWIEPH, latitude, longitude, system.swissCode, &cusps, &ascmc)
        if rc < 0 { throw SwissEphemerisError.calculationFailed(rc) }
        return HouseSet(
            system: system.rawValue,
            cusps: Array(cusps[1...12]),
            ascendant: ascmc[Int(SE_ASC)],
            mc: ascmc[Int(SE_MC)]
        )
    }

    private static func locateEphePath() -> String? {
        let fm = FileManager.default

        if let url = Bundle.main.url(forResource: "ephe", withExtension: nil),
           fm.fileExists(atPath: url.path) {
            return url.path
        }
        if let url = Bundle.main.url(forResource: "ephe", withExtension: nil, subdirectory: "Vendor/SwissEphemeris"),
           fm.fileExists(atPath: url.path) {
            return url.path
        }
        let local = URL(fileURLWithPath: "Vendor/SwissEphemeris/ephe", relativeTo: URL(fileURLWithPath: Bundle.main.bundlePath))
        if fm.fileExists(atPath: local.path) {
            return local.path
        }

        let sourcePath = URL(fileURLWithPath: #filePath)
            .deletingLastPathComponent()
            .appendingPathComponent("ephe")
        if fm.fileExists(atPath: sourcePath.path) {
            return sourcePath.path
        }
        return nil
    }
}

enum AstroMath {
    static func normalizedDegrees(_ value: Double) -> Double {
        var deg = value.truncatingRemainder(dividingBy: 360)
        if deg < 0 { deg += 360 }
        return deg
    }

    static func normalizedHours(_ value: Double) -> Double {
        var hours = value.truncatingRemainder(dividingBy: 24)
        if hours < 0 { hours += 24 }
        return hours
    }
}
