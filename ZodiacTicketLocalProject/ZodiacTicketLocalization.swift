import Foundation

public struct ZodiacTicketRuntimeTemplateItem: Codable, Equatable, Hashable {
    public let key: String
    public let text: String
    public let meta: [String: String]?

    public init(key: String, text: String, meta: [String: String]? = nil) {
        self.key = key
        self.text = text
        self.meta = meta
    }
}

public struct ZodiacTicketCopyTemplate: Codable, Equatable, Hashable {
    public let id: String
    public let sectionType: String
    public let domain: String
    public let intensityBucket: String
    public let headlineTemplate: String
    public let bodyTemplate: String
    public let timeHintTemplate: String?

    public init(
        id: String,
        sectionType: String,
        domain: String,
        intensityBucket: String,
        headlineTemplate: String,
        bodyTemplate: String,
        timeHintTemplate: String?
    ) {
        self.id = id
        self.sectionType = sectionType
        self.domain = domain
        self.intensityBucket = intensityBucket
        self.headlineTemplate = headlineTemplate
        self.bodyTemplate = bodyTemplate
        self.timeHintTemplate = timeHintTemplate
    }
}

public struct ZodiacTicketRuntimeDocument: Codable, Equatable {
    public let version: Int
    public let language: String
    public let strings: [String: String]
    public let templates: [String: [ZodiacTicketRuntimeTemplateItem]]
    public let astrologyCopyTemplates: [ZodiacTicketCopyTemplate]

    public init(
        version: Int,
        language: String,
        strings: [String: String],
        templates: [String: [ZodiacTicketRuntimeTemplateItem]],
        astrologyCopyTemplates: [ZodiacTicketCopyTemplate] = []
    ) {
        self.version = version
        self.language = language
        self.strings = strings
        self.templates = templates
        self.astrologyCopyTemplates = astrologyCopyTemplates
    }

    public static func empty(language: String) -> ZodiacTicketRuntimeDocument {
        ZodiacTicketRuntimeDocument(
            version: 1,
            language: language,
            strings: [:],
            templates: [:],
            astrologyCopyTemplates: []
        )
    }

    public func string(for key: String) -> String? {
        let value = strings[key]?.trimmingCharacters(in: .whitespacesAndNewlines)
        guard let value, !value.isEmpty else { return nil }
        return value
    }

    public func templateItems(for key: String) -> [ZodiacTicketRuntimeTemplateItem] {
        templates[key] ?? []
    }
}

public protocol ZodiacTicketTextProviding {
    func document(for language: ZodiacTicketLanguage) -> ZodiacTicketRuntimeDocument
}

private struct ZodiacTicketUnifiedDocument: Decodable {
    let version: Int?
    let languages: [String: ZodiacTicketUnifiedLanguagePayload]
}

private struct ZodiacTicketUnifiedLanguagePayload: Decodable {
    let language: String?
    let strings: [String: String]?
    let templates: [String: [ZodiacTicketRuntimeTemplateItem]]?
    let astrologyCopyTemplates: [ZodiacTicketCopyTemplate]?

    func document(version: Int, fallbackLanguage: String) -> ZodiacTicketRuntimeDocument {
        ZodiacTicketRuntimeDocument(
            version: version,
            language: language ?? fallbackLanguage,
            strings: strings ?? [:],
            templates: templates ?? [:],
            astrologyCopyTemplates: astrologyCopyTemplates ?? []
        )
    }
}

public final class ZodiacTicketBundleTextProvider: ZodiacTicketTextProviding {
    private let bundle: Bundle
    private let resourceDirectoryURL: URL?
    private let fileManager: FileManager
    private let decoder: JSONDecoder
    private var cache: [ZodiacTicketLanguage: ZodiacTicketRuntimeDocument]
    private let lock: NSLock

    public init(
        bundle: Bundle = .main,
        resourceDirectoryURL: URL? = nil,
        fileManager: FileManager = .default
    ) {
        self.bundle = bundle
        self.resourceDirectoryURL = resourceDirectoryURL
        self.fileManager = fileManager
        self.decoder = JSONDecoder()
        self.cache = [:]
        self.lock = NSLock()
    }

    public func document(for language: ZodiacTicketLanguage) -> ZodiacTicketRuntimeDocument {
        lock.lock()
        if let cached = cache[language] {
            lock.unlock()
            return cached
        }
        lock.unlock()

        let loaded = loadUnified(language: language)
            ?? loadRuntime(language: language)
            ?? (language == .en ? nil : loadRuntime(language: .en))
            ?? .empty(language: language.rawValue)

        lock.lock()
        cache[language] = loaded
        lock.unlock()
        return loaded
    }

    public func resetCache() {
        lock.lock()
        cache.removeAll()
        lock.unlock()
    }

    private func loadUnified(language: ZodiacTicketLanguage) -> ZodiacTicketRuntimeDocument? {
        for url in candidateURLs(fileName: "zodiac_ticket_text") {
            guard let data = try? Data(contentsOf: url),
                  let decoded = try? decoder.decode(ZodiacTicketUnifiedDocument.self, from: data) else {
                continue
            }

            for key in language.unifiedKeys {
                if let payload = decoded.languages[key] {
                    return payload.document(version: decoded.version ?? 1, fallbackLanguage: language.rawValue)
                }
            }
        }
        return nil
    }

    private func loadRuntime(language: ZodiacTicketLanguage) -> ZodiacTicketRuntimeDocument? {
        let fileName = "zodiac_ticket_runtime_\(language.runtimeResourceCode)"
        for url in candidateURLs(fileName: fileName) {
            guard let data = try? Data(contentsOf: url),
                  let decoded = try? decoder.decode(ZodiacTicketRuntimeDocument.self, from: data) else {
                continue
            }
            return decoded
        }
        return nil
    }

    private func candidateURLs(fileName: String) -> [URL] {
        let fileWithExtension = "\(fileName).json"
        var candidates: [URL] = []

        if let resourceDirectoryURL {
            candidates.append(resourceDirectoryURL.appendingPathComponent(fileWithExtension))
            candidates.append(resourceDirectoryURL.appendingPathComponent("Zodiac").appendingPathComponent(fileWithExtension))
            candidates.append(resourceDirectoryURL.appendingPathComponent("Resources").appendingPathComponent("Zodiac").appendingPathComponent(fileWithExtension))
        }

        let cwd = URL(fileURLWithPath: fileManager.currentDirectoryPath, isDirectory: true)
        candidates.append(cwd.appendingPathComponent(fileWithExtension))
        candidates.append(cwd.appendingPathComponent("Zodiac").appendingPathComponent(fileWithExtension))
        candidates.append(cwd.appendingPathComponent("Resources").appendingPathComponent("Zodiac").appendingPathComponent(fileWithExtension))
        candidates.append(cwd.appendingPathComponent("fDivinationApp").appendingPathComponent("Resources").appendingPathComponent("Zodiac").appendingPathComponent(fileWithExtension))

        if let resourceURL = bundle.resourceURL {
            candidates.append(resourceURL.appendingPathComponent(fileWithExtension))
            candidates.append(resourceURL.appendingPathComponent("Zodiac").appendingPathComponent(fileWithExtension))
            candidates.append(resourceURL.appendingPathComponent("Resources").appendingPathComponent("Zodiac").appendingPathComponent(fileWithExtension))
        }

        candidates.append(contentsOf: [
            bundle.url(forResource: fileName, withExtension: "json"),
            bundle.url(forResource: fileName, withExtension: "json", subdirectory: "Zodiac"),
            bundle.url(forResource: fileName, withExtension: "json", subdirectory: "Resources/Zodiac")
        ].compactMap { $0 })

        var seen = Set<String>()
        return candidates.filter { url in
            guard fileManager.fileExists(atPath: url.path), !seen.contains(url.path) else {
                return false
            }
            seen.insert(url.path)
            return true
        }
    }
}

public struct ZodiacTicketResolvedSegment: Codable, Equatable, Hashable {
    public let key: String
    public let tokens: [String: String]
    public let text: String

    public init(key: String, tokens: [String: String], text: String) {
        self.key = key
        self.tokens = tokens
        self.text = text
    }
}

public struct ZodiacTicketResolvedText: Codable, Equatable, Hashable {
    public let segments: [ZodiacTicketResolvedSegment]
    public let text: String

    public init(segments: [ZodiacTicketResolvedSegment], text: String) {
        self.segments = segments
        self.text = text
    }
}

public enum ZodiacTicketRuntimeRenderer {
    public static func renderNamed(
        _ item: ZodiacTicketRuntimeTemplateItem,
        tokens: [String: String]
    ) -> ZodiacTicketResolvedText {
        var rendered = item.text
        for (key, value) in tokens.sorted(by: { $0.key < $1.key }) {
            rendered = rendered.replacingOccurrences(of: "{\(key)}", with: value)
        }
        let segment = ZodiacTicketResolvedSegment(key: item.key, tokens: tokens, text: rendered)
        return ZodiacTicketResolvedText(segments: [segment], text: rendered)
    }

    public static func renderFormat(
        _ item: ZodiacTicketRuntimeTemplateItem,
        tokens: [String: String],
        orderedKeys: [String],
        locale: Locale
    ) -> ZodiacTicketResolvedText {
        let arguments = orderedKeys.map { tokens[$0] ?? "" }
        let rendered = String(format: item.text, locale: locale, arguments: arguments.map { $0 as CVarArg })
        let segment = ZodiacTicketResolvedSegment(key: item.key, tokens: tokens, text: rendered)
        return ZodiacTicketResolvedText(segments: [segment], text: rendered)
    }

    public static func renderLiteral(
        key: String,
        text: String,
        tokens: [String: String] = [:]
    ) -> ZodiacTicketResolvedText {
        ZodiacTicketResolvedText(
            segments: [ZodiacTicketResolvedSegment(key: key, tokens: tokens, text: text)],
            text: text
        )
    }
}

public enum ZodiacTicketLocalizationKey: String, CaseIterable {
    case title = "divination.zodiac.short"
    case serialNumber = "horoscope.ticket.label.serial_no"
    case inspiration = "horoscope.ticket.label.inspiration"
    case good = "horoscope.ticket.label.do"
    case caution = "horoscope.ticket.label.caution"
    case overview = "horoscope.ticket.label.overview"
    case triggerFocus = "horoscope.ticket.label.trigger_focus"
    case triggerBadge = "horoscope.ticket.label.trigger_badge"
    case moon = "horoscope.ticket.identity.moon"
    case asc = "horoscope.ticket.identity.asc"
    case moreFormat = "horoscope.ticket.plus_n_format"
}

public struct ZodiacTicketLocalizer {
    public let provider: ZodiacTicketTextProviding

    public init(provider: ZodiacTicketTextProviding = ZodiacTicketBundleTextProvider()) {
        self.provider = provider
    }

    public func text(_ key: String, language: ZodiacTicketLanguage) -> String {
        provider.document(for: language).string(for: key) ?? ""
    }

    public func labels(language: ZodiacTicketLanguage) -> ZodiacTicketLabels {
        ZodiacTicketLabels(
            title: text(ZodiacTicketLocalizationKey.title.rawValue, language: language),
            serialNumber: text(ZodiacTicketLocalizationKey.serialNumber.rawValue, language: language),
            inspiration: text(ZodiacTicketLocalizationKey.inspiration.rawValue, language: language),
            good: text(ZodiacTicketLocalizationKey.good.rawValue, language: language),
            caution: text(ZodiacTicketLocalizationKey.caution.rawValue, language: language),
            overview: text(ZodiacTicketLocalizationKey.overview.rawValue, language: language),
            triggerFocus: text(ZodiacTicketLocalizationKey.triggerFocus.rawValue, language: language),
            triggerBadge: text(ZodiacTicketLocalizationKey.triggerBadge.rawValue, language: language),
            moon: text(ZodiacTicketLocalizationKey.moon.rawValue, language: language),
            asc: text(ZodiacTicketLocalizationKey.asc.rawValue, language: language),
            moreFormat: text(ZodiacTicketLocalizationKey.moreFormat.rawValue, language: language)
        )
    }

    public func renderTemplate(
        key: String,
        language: ZodiacTicketLanguage,
        tokens: [String: String],
        selectionSeed: String
    ) -> ZodiacTicketResolvedText? {
        let items = provider.document(for: language).templateItems(for: key)
        guard !items.isEmpty else { return nil }
        let index = Int(ZodiacTicketUtilities.stableHash64(selectionSeed) % UInt64(items.count))
        return ZodiacTicketRuntimeRenderer.renderNamed(items[index], tokens: tokens)
    }
}

