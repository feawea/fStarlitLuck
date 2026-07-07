import Foundation

final class SwissEphAdapterImpl: SwissEphAdapter {
    private let service = SwissEphemerisService.shared

    func calcPlanet(jdUT: Double, planet: Planet) throws -> (lon: Double, speed: Double) {
        try service.calcPlanet(jdUT: jdUT, planet: planet)
    }

    func calcHouses(
        jdUT: Double,
        lat: Double,
        lon: Double,
        system: HouseSystem
    ) throws -> (cusps: [Double], asc: Double, mc: Double, armc: Double, trueObliquity: Double) {
        try service.calcHouseDetails(jdUT: jdUT, latitude: lat, longitude: lon, system: system)
    }
}
