import Foundation

public enum ZodiacTicketUtilities {
    public static func stableHash64(_ input: String) -> UInt64 {
        let prime: UInt64 = 1_099_511_628_211
        var hash: UInt64 = 14_695_981_039_346_656_037
        for byte in input.utf8 {
            hash ^= UInt64(byte)
            hash &*= prime
        }
        return hash
    }

    public static func normalizedInlineText(_ raw: String) -> String {
        raw.replacingOccurrences(of: "\n", with: " ")
            .replacingOccurrences(of: "\r", with: " ")
            .replacingOccurrences(of: "\u{2028}", with: " ")
            .replacingOccurrences(of: "\u{2029}", with: " ")
            .replacingOccurrences(of: "\u{200B}", with: "")
            .replacingOccurrences(of: "\u{200C}", with: "")
            .replacingOccurrences(of: "\u{200D}", with: "")
            .replacingOccurrences(of: "\u{2060}", with: "")
            .replacingOccurrences(of: "\u{00AD}", with: "")
            .replacingOccurrences(of: #"\s+"#, with: " ", options: .regularExpression)
            .trimmingCharacters(in: .whitespacesAndNewlines)
    }

    public static func compactDatePart(from storageKey: String) -> String? {
        let prefix = "astro_ticket_torn_"
        guard let range = storageKey.range(of: prefix) else { return nil }
        let suffix = storageKey[range.upperBound...]
        let digits = suffix.prefix(while: { $0.isNumber })
        guard digits.count >= 8 else { return nil }
        return String(digits.prefix(8))
    }

    public static func serialNumber(from storageKey: String) -> String {
        guard let datePart = compactDatePart(from: storageKey), datePart.count == 8 else {
            return ""
        }
        return "0\(datePart.dropFirst(2))"
    }

    public static func storageKey(dateKey: String, houseSystem: ZodiacTicketHouseSystem) -> String {
        let datePart = dateKey.replacingOccurrences(of: "-", with: "")
        return "astro_ticket_torn_\(datePart)_\(houseSystem.rawValue)"
    }
}

