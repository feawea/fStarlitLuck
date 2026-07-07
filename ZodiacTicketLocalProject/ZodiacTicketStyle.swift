import SwiftUI

public struct ZodiacTicketStyle {
    public let accent: Color
    public let paper: Color
    public let paperShadow: Color
    public let primaryText: Color
    public let bodyText: Color
    public let secondaryText: Color
    public let mutedText: Color
    public let separator: Color
    public let cornerRadius: CGFloat

    public init(
        accent: Color = Color(red: 0.23, green: 0.36, blue: 0.64),
        paper: Color = Color(red: 0.98, green: 0.965, blue: 0.92),
        paperShadow: Color = Color.black.opacity(0.16),
        primaryText: Color = Color(red: 0.13, green: 0.18, blue: 0.28),
        bodyText: Color = Color(red: 0.18, green: 0.23, blue: 0.34),
        secondaryText: Color = Color(red: 0.38, green: 0.45, blue: 0.58),
        mutedText: Color = Color(red: 0.53, green: 0.60, blue: 0.70),
        separator: Color = Color(red: 0.23, green: 0.36, blue: 0.64).opacity(0.16),
        cornerRadius: CGFloat = 18
    ) {
        self.accent = accent
        self.paper = paper
        self.paperShadow = paperShadow
        self.primaryText = primaryText
        self.bodyText = bodyText
        self.secondaryText = secondaryText
        self.mutedText = mutedText
        self.separator = separator
        self.cornerRadius = cornerRadius
    }
}

private struct ZodiacTicketDebossedTextModifier: ViewModifier {
    let enabled: Bool
    let baseColor: Color
    let depth: CGFloat
    let inkOpacity: Double

    func body(content: Content) -> some View {
        if enabled {
            content
                .foregroundStyle(baseColor.opacity(inkOpacity))
                .shadow(color: .white.opacity(0.55), radius: 0, x: 0, y: -depth * 0.38)
                .shadow(color: .black.opacity(0.10), radius: 0, x: 0, y: depth * 0.42)
        } else {
            content.foregroundStyle(baseColor)
        }
    }
}

extension View {
    func zodiacTicketDebossedText(
        enabled: Bool,
        baseColor: Color,
        depth: CGFloat = 1.4,
        inkOpacity: Double = 0.78
    ) -> some View {
        modifier(
            ZodiacTicketDebossedTextModifier(
                enabled: enabled,
                baseColor: baseColor,
                depth: depth,
                inkOpacity: inkOpacity
            )
        )
    }
}
