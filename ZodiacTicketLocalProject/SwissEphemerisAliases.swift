import Foundation

typealias Planet = ZodiacTicketPlanet

extension ZodiacTicketPlanet {
    var swissCode: Int32 {
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

