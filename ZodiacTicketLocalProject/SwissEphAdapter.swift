import Foundation

protocol SwissEphAdapter {
    func calcPlanet(jdUT: Double, planet: Planet) throws -> (lon: Double, speed: Double)
    func calcHouses(
        jdUT: Double,
        lat: Double,
        lon: Double,
        system: HouseSystem
    ) throws -> (cusps: [Double], asc: Double, mc: Double, armc: Double, trueObliquity: Double)
}
