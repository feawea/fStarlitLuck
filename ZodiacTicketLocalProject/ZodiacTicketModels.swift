import Foundation

public enum ZodiacTicketPlanet: String, CaseIterable, Codable, Hashable {
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

    public var sortIndex: Int {
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

public enum ZodiacTicketAspect: String, CaseIterable, Codable, Hashable {
    case conjunction
    case sextile
    case square
    case trine
    case opposition

    public var angle: Double {
        switch self {
        case .conjunction: return 0
        case .sextile: return 60
        case .square: return 90
        case .trine: return 120
        case .opposition: return 180
        }
    }
}

public enum ZodiacTicketDomain: String, CaseIterable, Codable, Hashable {
    case selfDomain = "self"
    case money
    case communication
    case home
    case loveCreativity = "love_creativity"
    case healthWorkflow = "health_workflow"
    case relationship
    case riskShared = "risk_shared"
    case travelLearning = "travel_learning"
    case career
    case socialFuture = "social_future"
    case restInner = "rest_inner"
}

public enum ZodiacTicketHouseSystem: String, CaseIterable, Codable, Hashable, Identifiable {
    case placidus = "P"
    case koch = "K"
    case equal = "E"

    public var id: String { rawValue }
}

public enum ZodiacTicketTarget: Codable, Hashable {
    case planet(ZodiacTicketPlanet)
    case asc
    case mc

    private enum CodingKeys: String, CodingKey {
        case kind
        case planet
    }

    private enum Kind: String, Codable {
        case planet
        case asc
        case mc
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        let kind = try container.decode(Kind.self, forKey: .kind)
        switch kind {
        case .planet:
            self = .planet(try container.decode(ZodiacTicketPlanet.self, forKey: .planet))
        case .asc:
            self = .asc
        case .mc:
            self = .mc
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        switch self {
        case .planet(let planet):
            try container.encode(Kind.planet, forKey: .kind)
            try container.encode(planet, forKey: .planet)
        case .asc:
            try container.encode(Kind.asc, forKey: .kind)
        case .mc:
            try container.encode(Kind.mc, forKey: .kind)
        }
    }
}

public struct ZodiacTicketTransitEvent: Codable, Hashable {
    public let hour: Int
    public let transitPlanet: ZodiacTicketPlanet
    public let target: ZodiacTicketTarget
    public let aspect: ZodiacTicketAspect
    public let orbDegrees: Double
    public let applying: Bool
    public let transitHouse: Int
    public let targetHouse: Int
    public let score: Double
    public let direction: Double

    public init(
        hour: Int,
        transitPlanet: ZodiacTicketPlanet,
        target: ZodiacTicketTarget,
        aspect: ZodiacTicketAspect,
        orbDegrees: Double,
        applying: Bool,
        transitHouse: Int,
        targetHouse: Int,
        score: Double,
        direction: Double
    ) {
        self.hour = hour
        self.transitPlanet = transitPlanet
        self.target = target
        self.aspect = aspect
        self.orbDegrees = orbDegrees
        self.applying = applying
        self.transitHouse = transitHouse
        self.targetHouse = targetHouse
        self.score = score
        self.direction = direction
    }
}

public struct ZodiacTicketPeak: Codable, Hashable {
    public let domain: ZodiacTicketDomain
    public let peakLocalHour: Int
    public let peakScore: Double
    public let polarity: String
    public let intensity: Double
    public let drivers: [ZodiacTicketTransitEvent]

    public init(
        domain: ZodiacTicketDomain,
        peakLocalHour: Int,
        peakScore: Double,
        polarity: String,
        intensity: Double,
        drivers: [ZodiacTicketTransitEvent]
    ) {
        self.domain = domain
        self.peakLocalHour = peakLocalHour
        self.peakScore = peakScore
        self.polarity = polarity
        self.intensity = intensity
        self.drivers = drivers
    }
}

public struct ZodiacTicketExtras: Codable, Hashable {
    public let luckyWindow: String
    public let luckyPlanet: ZodiacTicketPlanet
    public let luckyNumber: Int
    public let luckyColorKey: String
    public let keywordKey: String

    public init(
        luckyWindow: String,
        luckyPlanet: ZodiacTicketPlanet,
        luckyNumber: Int,
        luckyColorKey: String,
        keywordKey: String
    ) {
        self.luckyWindow = luckyWindow
        self.luckyPlanet = luckyPlanet
        self.luckyNumber = luckyNumber
        self.luckyColorKey = luckyColorKey
        self.keywordKey = keywordKey
    }
}

public struct ZodiacTicketInput: Codable, Hashable {
    public let dateKey: String
    public let timeZoneIdentifier: String
    public let houseSystem: ZodiacTicketHouseSystem
    public let reminder: ZodiacTicketPeak
    public let goodLuck: [ZodiacTicketPeak]
    public let caution: [ZodiacTicketPeak]
    public let extras: ZodiacTicketExtras

    public init(
        dateKey: String,
        timeZoneIdentifier: String,
        houseSystem: ZodiacTicketHouseSystem,
        reminder: ZodiacTicketPeak,
        goodLuck: [ZodiacTicketPeak],
        caution: [ZodiacTicketPeak],
        extras: ZodiacTicketExtras
    ) {
        self.dateKey = dateKey
        self.timeZoneIdentifier = timeZoneIdentifier
        self.houseSystem = houseSystem
        self.reminder = reminder
        self.goodLuck = goodLuck
        self.caution = caution
        self.extras = extras
    }
}

public struct ZodiacTicketBigThreeItem: Codable, Hashable, Identifiable {
    public let id: String
    public let iconSystemName: String
    public let shortLabel: String
    public let signName: String
    public let oneLineHint: String

    public init(
        id: String,
        iconSystemName: String,
        shortLabel: String,
        signName: String,
        oneLineHint: String
    ) {
        self.id = id
        self.iconSystemName = iconSystemName
        self.shortLabel = shortLabel
        self.signName = signName
        self.oneLineHint = oneLineHint
    }
}

public struct ZodiacTicketBigThree: Codable, Hashable {
    public let sun: ZodiacTicketBigThreeItem
    public let moon: ZodiacTicketBigThreeItem
    public let asc: ZodiacTicketBigThreeItem

    public init(
        sun: ZodiacTicketBigThreeItem,
        moon: ZodiacTicketBigThreeItem,
        asc: ZodiacTicketBigThreeItem
    ) {
        self.sun = sun
        self.moon = moon
        self.asc = asc
    }
}

public struct ZodiacTicketCompactFields: Codable, Hashable {
    public let sunMoonAscLine: String
    public let inspirationShort: String
    public let luckyTitleShort: String
    public let luckyOneLineShort: String
    public let luckyMoreCount: Int
    public let cautionTitleShort: String
    public let cautionOneLineShort: String
    public let cautionMoreCount: Int
    public let focusSummaryShort: String
    public let curiosityQuestionShort: String
    public let plainReasonShort: String
    public let triggerHookLine: String
    public let revealLine1: String
    public let revealLine2: String?
    public let didTearTodayKey: String
    public let stampAssetName: String
    public let stampOpacity: Double
    public let stampRotationDegrees: Double

    public init(
        sunMoonAscLine: String,
        inspirationShort: String,
        luckyTitleShort: String,
        luckyOneLineShort: String,
        luckyMoreCount: Int,
        cautionTitleShort: String,
        cautionOneLineShort: String,
        cautionMoreCount: Int,
        focusSummaryShort: String,
        curiosityQuestionShort: String,
        plainReasonShort: String,
        triggerHookLine: String,
        revealLine1: String,
        revealLine2: String?,
        didTearTodayKey: String,
        stampAssetName: String,
        stampOpacity: Double,
        stampRotationDegrees: Double
    ) {
        self.sunMoonAscLine = sunMoonAscLine
        self.inspirationShort = inspirationShort
        self.luckyTitleShort = luckyTitleShort
        self.luckyOneLineShort = luckyOneLineShort
        self.luckyMoreCount = luckyMoreCount
        self.cautionTitleShort = cautionTitleShort
        self.cautionOneLineShort = cautionOneLineShort
        self.cautionMoreCount = cautionMoreCount
        self.focusSummaryShort = focusSummaryShort
        self.curiosityQuestionShort = curiosityQuestionShort
        self.plainReasonShort = plainReasonShort
        self.triggerHookLine = triggerHookLine
        self.revealLine1 = revealLine1
        self.revealLine2 = revealLine2
        self.didTearTodayKey = didTearTodayKey
        self.stampAssetName = stampAssetName
        self.stampOpacity = stampOpacity
        self.stampRotationDegrees = stampRotationDegrees
    }
}

public struct ZodiacTicketLabels: Codable, Hashable {
    public let title: String
    public let serialNumber: String
    public let inspiration: String
    public let good: String
    public let caution: String
    public let overview: String
    public let triggerFocus: String
    public let triggerBadge: String
    public let moon: String
    public let asc: String
    public let moreFormat: String

    public init(
        title: String,
        serialNumber: String,
        inspiration: String,
        good: String,
        caution: String,
        overview: String,
        triggerFocus: String,
        triggerBadge: String,
        moon: String,
        asc: String,
        moreFormat: String
    ) {
        self.title = title
        self.serialNumber = serialNumber
        self.inspiration = inspiration
        self.good = good
        self.caution = caution
        self.overview = overview
        self.triggerFocus = triggerFocus
        self.triggerBadge = triggerBadge
        self.moon = moon
        self.asc = asc
        self.moreFormat = moreFormat
    }
}

public struct ZodiacTicketViewModel: Codable, Hashable {
    public let language: ZodiacTicketLanguage
    public let displayName: String
    public let labels: ZodiacTicketLabels
    public let bigThree: ZodiacTicketBigThree?
    public let compactFields: ZodiacTicketCompactFields
    public let serialNumber: String
    public let tearCount: Int

    public init(
        language: ZodiacTicketLanguage,
        displayName: String,
        labels: ZodiacTicketLabels,
        bigThree: ZodiacTicketBigThree?,
        compactFields: ZodiacTicketCompactFields,
        serialNumber: String,
        tearCount: Int
    ) {
        self.language = language
        self.displayName = displayName
        self.labels = labels
        self.bigThree = bigThree
        self.compactFields = compactFields
        self.serialNumber = serialNumber
        self.tearCount = tearCount
    }
}

