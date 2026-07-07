import Foundation

public enum ZodiacTicketCompactTemplateKey: String, CaseIterable {
    case sunMoonAscLine = "ticket.compact.sun_moon_asc"
    case inspirationShort = "ticket.compact.inspiration"
    case luckyTitleShort = "ticket.compact.lucky_title"
    case luckyOneLineShort = "ticket.compact.lucky_line"
    case cautionTitleShort = "ticket.compact.caution_title"
    case cautionOneLineShort = "ticket.compact.caution_line"
    case focusSummaryShort = "ticket.compact.focus_summary"
    case curiosityQuestionShort = "ticket.compact.curiosity_question"
    case plainReasonShort = "ticket.compact.plain_reason"
    case triggerHookLine = "ticket.compact.trigger_hook"
    case revealLine1 = "ticket.compact.reveal_line_1"
    case revealLine2 = "ticket.compact.reveal_line_2"
}

public struct ZodiacTicketComposer {
    public let localizer: ZodiacTicketLocalizer

    public init(localizer: ZodiacTicketLocalizer = ZodiacTicketLocalizer()) {
        self.localizer = localizer
    }

    public func makeCompactFields(
        input: ZodiacTicketInput,
        bigThree: ZodiacTicketBigThree?,
        language: ZodiacTicketLanguage,
        userSeed: String
    ) -> ZodiacTicketCompactFields {
        let storageKey = ZodiacTicketUtilities.storageKey(dateKey: input.dateKey, houseSystem: input.houseSystem)
        let tokens = baseTokens(input: input, bigThree: bigThree, language: language, userSeed: userSeed)
        let datePart = input.dateKey.replacingOccurrences(of: "-", with: "")

        return ZodiacTicketCompactFields(
            sunMoonAscLine: rendered(.sunMoonAscLine, language: language, tokens: tokens, seed: "\(storageKey)|sunMoonAsc"),
            inspirationShort: rendered(.inspirationShort, language: language, tokens: tokens, seed: "\(storageKey)|inspiration"),
            luckyTitleShort: rendered(.luckyTitleShort, language: language, tokens: tokens, seed: "\(storageKey)|luckyTitle"),
            luckyOneLineShort: rendered(.luckyOneLineShort, language: language, tokens: tokens, seed: "\(storageKey)|luckyLine"),
            luckyMoreCount: max(0, input.goodLuck.count - 1),
            cautionTitleShort: rendered(.cautionTitleShort, language: language, tokens: tokens, seed: "\(storageKey)|cautionTitle"),
            cautionOneLineShort: rendered(.cautionOneLineShort, language: language, tokens: tokens, seed: "\(storageKey)|cautionLine"),
            cautionMoreCount: max(0, input.caution.count - 1),
            focusSummaryShort: rendered(.focusSummaryShort, language: language, tokens: tokens, seed: "\(storageKey)|focus"),
            curiosityQuestionShort: rendered(.curiosityQuestionShort, language: language, tokens: tokens, seed: "\(storageKey)|question"),
            plainReasonShort: rendered(.plainReasonShort, language: language, tokens: tokens, seed: "\(storageKey)|reason"),
            triggerHookLine: rendered(.triggerHookLine, language: language, tokens: tokens, seed: "\(storageKey)|hook"),
            revealLine1: rendered(.revealLine1, language: language, tokens: tokens, seed: "\(storageKey)|reveal1"),
            revealLine2: optionalRendered(.revealLine2, language: language, tokens: tokens, seed: "\(storageKey)|reveal2"),
            didTearTodayKey: storageKey,
            stampAssetName: stampAssetName(datePart: datePart, houseSystem: input.houseSystem, userSeed: userSeed),
            stampOpacity: stampOpacity(datePart: datePart, houseSystem: input.houseSystem, userSeed: userSeed),
            stampRotationDegrees: stampRotation(datePart: datePart, houseSystem: input.houseSystem, userSeed: userSeed)
        )
    }

    public func makeViewModel(
        compactFields: ZodiacTicketCompactFields,
        bigThree: ZodiacTicketBigThree?,
        displayName: String,
        language: ZodiacTicketLanguage,
        tearCount: Int
    ) -> ZodiacTicketViewModel {
        ZodiacTicketViewModel(
            language: language,
            displayName: displayName.trimmingCharacters(in: .whitespacesAndNewlines),
            labels: localizer.labels(language: language),
            bigThree: bigThree,
            compactFields: compactFields,
            serialNumber: ZodiacTicketUtilities.serialNumber(from: compactFields.didTearTodayKey),
            tearCount: max(0, tearCount)
        )
    }

    public func makeViewModel(
        input: ZodiacTicketInput,
        bigThree: ZodiacTicketBigThree?,
        displayName: String,
        language: ZodiacTicketLanguage,
        userSeed: String,
        tearCount: Int
    ) -> ZodiacTicketViewModel {
        let fields = makeCompactFields(
            input: input,
            bigThree: bigThree,
            language: language,
            userSeed: userSeed
        )
        return makeViewModel(
            compactFields: fields,
            bigThree: bigThree,
            displayName: displayName,
            language: language,
            tearCount: tearCount
        )
    }

    private func rendered(
        _ key: ZodiacTicketCompactTemplateKey,
        language: ZodiacTicketLanguage,
        tokens: [String: String],
        seed: String
    ) -> String {
        localizer.renderTemplate(
            key: key.rawValue,
            language: language,
            tokens: tokens,
            selectionSeed: seed
        )?.text ?? ""
    }

    private func optionalRendered(
        _ key: ZodiacTicketCompactTemplateKey,
        language: ZodiacTicketLanguage,
        tokens: [String: String],
        seed: String
    ) -> String? {
        let value = rendered(key, language: language, tokens: tokens, seed: seed)
        return value.isEmpty ? nil : value
    }

    private func baseTokens(
        input: ZodiacTicketInput,
        bigThree: ZodiacTicketBigThree?,
        language: ZodiacTicketLanguage,
        userSeed: String
    ) -> [String: String] {
        var tokens: [String: String] = [
            "date_key": input.dateKey,
            "timezone": input.timeZoneIdentifier,
            "house_system": input.houseSystem.rawValue,
            "language": language.rawValue,
            "user_seed": userSeed,
            "reminder_domain": input.reminder.domain.rawValue,
            "reminder_hour": String(input.reminder.peakLocalHour),
            "reminder_intensity": String(format: "%.2f", input.reminder.intensity),
            "good_domain": input.goodLuck.first?.domain.rawValue ?? "",
            "good_hour": input.goodLuck.first.map { String($0.peakLocalHour) } ?? "",
            "caution_domain": input.caution.first?.domain.rawValue ?? "",
            "caution_hour": input.caution.first.map { String($0.peakLocalHour) } ?? "",
            "lucky_window": input.extras.luckyWindow,
            "lucky_planet": input.extras.luckyPlanet.rawValue,
            "lucky_number": String(input.extras.luckyNumber),
            "lucky_color_key": input.extras.luckyColorKey,
            "keyword_key": input.extras.keywordKey
        ]

        if let bigThree {
            tokens["sun_label"] = bigThree.sun.shortLabel
            tokens["sun_sign"] = bigThree.sun.signName
            tokens["moon_label"] = bigThree.moon.shortLabel
            tokens["moon_sign"] = bigThree.moon.signName
            tokens["asc_label"] = bigThree.asc.shortLabel
            tokens["asc_sign"] = bigThree.asc.signName
        }

        if let driver = strongestDriver(from: input.reminder.drivers) {
            tokens["trigger_planet"] = driver.transitPlanet.rawValue
            tokens["trigger_aspect"] = driver.aspect.rawValue
            tokens["trigger_target"] = targetToken(driver.target)
            tokens["trigger_house"] = String(driver.transitHouse)
            tokens["trigger_orb"] = String(format: "%.2f", driver.orbDegrees)
        }

        return tokens
    }

    private func strongestDriver(from drivers: [ZodiacTicketTransitEvent]) -> ZodiacTicketTransitEvent? {
        drivers.sorted {
            if $0.score == $1.score {
                return $0.hour < $1.hour
            }
            return $0.score > $1.score
        }.first
    }

    private func targetToken(_ target: ZodiacTicketTarget) -> String {
        switch target {
        case .planet(let planet): return planet.rawValue
        case .asc: return "asc"
        case .mc: return "mc"
        }
    }

    private func stampAssetName(datePart: String, houseSystem: ZodiacTicketHouseSystem, userSeed: String) -> String {
        let seed = "\(datePart)|\(houseSystem.rawValue)|\(userSeed)|stamp"
        let index = Int(ZodiacTicketUtilities.stableHash64(seed) % 5) + 1
        switch index {
        case 1: return "stamp_1_red"
        case 2: return "stamp_2_green"
        case 3: return "stamp_3_blue"
        case 4: return "stamp_4_brown"
        default: return "stamp_5_oval"
        }
    }

    private func stampOpacity(datePart: String, houseSystem: ZodiacTicketHouseSystem, userSeed: String) -> Double {
        let seed = "\(datePart)|\(houseSystem.rawValue)|\(userSeed)|stampOpacity"
        return 0.12 + Double(ZodiacTicketUtilities.stableHash64(seed) % 18) / 100
    }

    private func stampRotation(datePart: String, houseSystem: ZodiacTicketHouseSystem, userSeed: String) -> Double {
        let seed = "\(datePart)|\(houseSystem.rawValue)|\(userSeed)|stampRotation"
        return Double(Int(ZodiacTicketUtilities.stableHash64(seed) % 17)) - 8
    }
}

