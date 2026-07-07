import Foundation

public enum ZodiacTicketLanguage: String, CaseIterable, Codable, Hashable, Identifiable {
    case zhHans = "zh-Hans"
    case en = "en"
    case zhHant = "zh-Hant"
    case ja = "ja"
    case de = "de"
    case es = "es"
    case pt = "pt"
    case fr = "fr"
    case ko = "ko"
    case ru = "ru"
    case it = "it"

    public var id: String { rawValue }

    public var runtimeResourceCode: String {
        switch self {
        case .zhHans: return "zh_hans"
        case .en: return "en_us"
        case .zhHant: return "zh_hant"
        case .ja: return "ja_jp"
        case .de: return "de_de"
        case .es: return "es_es"
        case .pt: return "pt_pt"
        case .fr: return "fr_fr"
        case .ko: return "ko_kr"
        case .ru: return "ru_ru"
        case .it: return "it_it"
        }
    }

    public var unifiedKeys: [String] {
        [runtimeResourceCode, rawValue, rawValue.replacingOccurrences(of: "-", with: "_")]
    }

    public static func resolve(from locale: Locale) -> ZodiacTicketLanguage {
        let normalized = locale.identifier
            .replacingOccurrences(of: "_", with: "-")
            .lowercased()

        if normalized.contains("zh-hant") || normalized.hasPrefix("zh-tw") || normalized.hasPrefix("zh-hk") || normalized.hasPrefix("zh-mo") {
            return .zhHant
        }
        if normalized == "zh" || normalized.hasPrefix("zh-") {
            return .zhHans
        }
        if normalized == "ja" || normalized.hasPrefix("ja-") { return .ja }
        if normalized == "de" || normalized.hasPrefix("de-") { return .de }
        if normalized == "es" || normalized.hasPrefix("es-") { return .es }
        if normalized == "pt" || normalized.hasPrefix("pt-") { return .pt }
        if normalized == "fr" || normalized.hasPrefix("fr-") { return .fr }
        if normalized == "ko" || normalized.hasPrefix("ko-") { return .ko }
        if normalized == "ru" || normalized.hasPrefix("ru-") { return .ru }
        if normalized == "it" || normalized.hasPrefix("it-") { return .it }
        return .en
    }
}

