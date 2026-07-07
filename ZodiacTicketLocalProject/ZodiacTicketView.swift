import SwiftUI

public struct ZodiacTicketView: View {
    public let model: ZodiacTicketViewModel
    public let style: ZodiacTicketStyle
    public let usesDebossedText: Bool

    public init(
        model: ZodiacTicketViewModel,
        style: ZodiacTicketStyle = ZodiacTicketStyle(),
        usesDebossedText: Bool = true
    ) {
        self.model = model
        self.style = style
        self.usesDebossedText = usesDebossedText
    }

    public var body: some View {
        VStack(alignment: .leading, spacing: 0) {
            upperSection
            divider
                .padding(.vertical, 10)
            lowerSection
        }
        .padding(.horizontal, 18)
        .padding(.vertical, 18)
        .background(ticketBackground)
        .clipShape(RoundedRectangle(cornerRadius: style.cornerRadius, style: .continuous))
        .shadow(color: style.paperShadow, radius: 16, x: 0, y: 8)
        .accessibilityElement(children: .contain)
    }

    private var upperSection: some View {
        VStack(alignment: .leading, spacing: 12) {
            HStack(alignment: .top, spacing: 12) {
                VStack(alignment: .leading, spacing: 8) {
                    optionalText(model.labels.title)
                        .font(.system(size: languageUsesCJK ? 24 : 22, weight: .black))
                        .zodiacTicketDebossedText(
                            enabled: usesDebossedText,
                            baseColor: style.primaryText,
                            depth: 1.8,
                            inkOpacity: 0.90
                        )
                        .lineLimit(languageUsesCJK ? 1 : 2)
                        .minimumScaleFactor(0.74)

                    Rectangle()
                        .fill(style.accent.opacity(0.22))
                        .frame(width: 116, height: 2)
                }
                .background(alignment: .topLeading) {
                    Text(String(model.tearCount))
                        .font(.system(size: 76, weight: .bold, design: .rounded))
                        .foregroundStyle(style.accent.opacity(0.10))
                        .offset(x: 22, y: -16)
                        .accessibilityHidden(true)
                }

                Spacer(minLength: 8)

                serialStamp
            }

            identityBlock
        }
    }

    private var identityBlock: some View {
        VStack(alignment: .leading, spacing: 6) {
            if let bigThree = model.bigThree {
                HStack(alignment: .firstTextBaseline, spacing: 6) {
                    optionalText(model.displayName)
                        .font(.system(size: 15, weight: .bold))
                    optionalText(bigThree.sun.signName)
                        .font(.system(size: 15, weight: .semibold))
                }
                .zodiacTicketDebossedText(
                    enabled: usesDebossedText,
                    baseColor: style.bodyText,
                    depth: 1.35,
                    inkOpacity: 0.82
                )

                HStack(alignment: .firstTextBaseline, spacing: 6) {
                    optionalText(model.labels.moon)
                        .font(.system(size: 12, weight: .semibold))
                    optionalText(bigThree.moon.signName)
                        .font(.system(size: 12, weight: .regular))
                    optionalText(model.labels.asc)
                        .font(.system(size: 12, weight: .semibold))
                    optionalText(bigThree.asc.signName)
                        .font(.system(size: 12, weight: .regular))
                }
                .zodiacTicketDebossedText(
                    enabled: usesDebossedText,
                    baseColor: style.secondaryText,
                    depth: 1.15,
                    inkOpacity: 0.72
                )
            } else {
                optionalText(model.compactFields.sunMoonAscLine)
                    .font(.system(size: 13, weight: .semibold))
                    .zodiacTicketDebossedText(
                        enabled: usesDebossedText,
                        baseColor: style.bodyText,
                        depth: 1.25,
                        inkOpacity: 0.78
                    )
            }
        }
        .padding(.vertical, 8)
        .frame(maxWidth: .infinity, alignment: .leading)
        .overlay(alignment: .top) {
            Rectangle()
                .fill(style.separator)
                .frame(height: 1)
        }
        .overlay(alignment: .bottom) {
            Rectangle()
                .fill(style.separator)
                .frame(height: 1)
        }
    }

    private var lowerSection: some View {
        VStack(alignment: .leading, spacing: 12) {
            headlineBlock

            row(
                iconSystemName: "sparkles",
                title: model.labels.inspiration,
                text: model.compactFields.inspirationShort,
                moreCount: nil
            )

            divider

            row(
                iconSystemName: "leaf.fill",
                title: model.compactFields.luckyTitleShort.isEmpty ? model.labels.good : model.compactFields.luckyTitleShort,
                text: model.compactFields.luckyOneLineShort,
                moreCount: model.compactFields.luckyMoreCount
            )

            divider

            row(
                iconSystemName: "exclamationmark.triangle.fill",
                title: model.compactFields.cautionTitleShort.isEmpty ? model.labels.caution : model.compactFields.cautionTitleShort,
                text: model.compactFields.cautionOneLineShort,
                moreCount: model.compactFields.cautionMoreCount
            )

            triggerBlock
        }
    }

    private var headlineBlock: some View {
        VStack(alignment: .leading, spacing: 7) {
            optionalText(model.labels.overview)
                .font(.system(size: 12, weight: .semibold))
                .foregroundStyle(style.secondaryText)

            optionalText(model.compactFields.focusSummaryShort)
                .font(.system(size: languageUsesCJK ? 22 : 20, weight: .bold))
                .zodiacTicketDebossedText(
                    enabled: usesDebossedText,
                    baseColor: style.primaryText,
                    depth: 1.6,
                    inkOpacity: 0.88
                )
                .lineLimit(3)
                .minimumScaleFactor(0.82)
                .fixedSize(horizontal: false, vertical: true)
        }
        .padding(.horizontal, 8)
        .padding(.vertical, 10)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(style.accent.opacity(0.07), in: RoundedRectangle(cornerRadius: 10, style: .continuous))
    }

    private var triggerBlock: some View {
        VStack(alignment: .leading, spacing: 7) {
            HStack(spacing: 8) {
                optionalText(model.labels.triggerFocus)
                    .font(.system(size: 12, weight: .semibold))
                    .foregroundStyle(style.secondaryText)

                Spacer(minLength: 0)

                optionalText(model.labels.triggerBadge)
                    .font(.system(size: 11, weight: .semibold))
                    .foregroundStyle(style.secondaryText)
                    .padding(.horizontal, 7)
                    .padding(.vertical, 2)
                    .background(style.mutedText.opacity(0.12), in: Capsule())
            }

            optionalText(model.compactFields.revealLine1)
                .font(.system(size: 15, weight: .regular))
                .zodiacTicketDebossedText(
                    enabled: usesDebossedText,
                    baseColor: style.bodyText,
                    depth: 1.35,
                    inkOpacity: 0.74
                )
                .lineLimit(nil)
                .fixedSize(horizontal: false, vertical: true)

            if let revealLine2 = model.compactFields.revealLine2, !revealLine2.isEmpty {
                optionalText(revealLine2)
                    .font(.system(size: 14, weight: .regular))
                    .zodiacTicketDebossedText(
                        enabled: usesDebossedText,
                        baseColor: style.secondaryText,
                        depth: 1.2,
                        inkOpacity: 0.64
                    )
                    .lineLimit(nil)
                    .fixedSize(horizontal: false, vertical: true)
            }
        }
        .padding(.horizontal, 10)
        .padding(.vertical, 10)
        .background(Color.white.opacity(0.20), in: RoundedRectangle(cornerRadius: 10, style: .continuous))
        .overlay(
            RoundedRectangle(cornerRadius: 10, style: .continuous)
                .stroke(style.accent.opacity(0.18), lineWidth: 1)
        )
    }

    private func row(
        iconSystemName: String,
        title: String,
        text: String,
        moreCount: Int?
    ) -> some View {
        HStack(alignment: .top, spacing: 10) {
            Image(systemName: iconSystemName)
                .font(.system(size: 14, weight: .semibold))
                .foregroundStyle(style.accent)
                .frame(width: 18)
                .padding(.top, 2)
                .accessibilityHidden(true)

            VStack(alignment: .leading, spacing: 3) {
                optionalText(title)
                    .font(.system(size: 12, weight: .semibold))
                    .foregroundStyle(style.secondaryText)

                optionalText(text)
                    .font(.system(size: 15, weight: .medium))
                    .zodiacTicketDebossedText(
                        enabled: usesDebossedText,
                        baseColor: style.bodyText,
                        depth: 1.35,
                        inkOpacity: 0.76
                    )
                    .lineLimit(2)
                    .minimumScaleFactor(0.88)
                    .fixedSize(horizontal: false, vertical: true)
            }
            .frame(maxWidth: .infinity, alignment: .leading)

            if let moreCount, moreCount > 0, let moreText = moreChipText(count: moreCount) {
                Text(moreText)
                    .font(.system(size: 12, weight: .semibold))
                    .foregroundStyle(style.secondaryText)
                    .padding(.horizontal, 7)
                    .padding(.vertical, 3)
                    .background(style.mutedText.opacity(0.12), in: Capsule())
            }
        }
        .frame(maxWidth: .infinity, alignment: .topLeading)
    }

    private var serialStamp: some View {
        VStack(alignment: .trailing, spacing: 3) {
            optionalText(model.labels.serialNumber)
                .font(.system(size: 10, weight: .bold))
                .foregroundStyle(style.secondaryText)
                .kerning(1.0)

            optionalText(model.serialNumber)
                .font(.system(size: 14, weight: .semibold, design: .rounded).monospacedDigit())
                .foregroundStyle(style.bodyText)
                .kerning(0.8)

            barcode(serialNumber: model.serialNumber)
                .frame(width: 76, height: 7)
                .padding(.top, 2)
        }
        .fixedSize()
    }

    private func barcode(serialNumber: String) -> some View {
        GeometryReader { proxy in
            let segments = barcodeSegments(serialNumber: serialNumber)
            let totalUnits = max(1, segments.reduce(0) { $0 + $1.units })
            let unit = proxy.size.width / CGFloat(totalUnits)
            let barInset = max(0, unit * 0.08)

            Path { path in
                var x: CGFloat = 0
                for segment in segments {
                    let segmentWidth = CGFloat(segment.units) * unit
                    if segment.isBar {
                        path.addRect(
                            CGRect(
                                x: x + barInset,
                                y: 0,
                                width: max(0.5, segmentWidth - barInset * 2),
                                height: proxy.size.height
                            )
                        )
                    }
                    x += segmentWidth
                }
            }
            .fill(style.bodyText)
        }
        .accessibilityHidden(true)
    }

    private func barcodeSegments(serialNumber: String) -> [(isBar: Bool, units: Int)] {
        let digits = Array(serialNumber.filter(\.isNumber))
        var segments: [(isBar: Bool, units: Int)] = []
        var isBar = true

        func append(_ units: Int) {
            segments.append((isBar: isBar, units: max(1, units)))
            isBar.toggle()
        }

        [2, 1, 1, 1, 2, 1, 1].forEach(append)

        if digits.isEmpty {
            [1, 1, 2, 1, 1, 2, 1, 1].forEach(append)
        } else {
            for (index, character) in digits.enumerated() {
                let value = Int(String(character)) ?? 0
                [
                    1 + ((value + index) % 3),
                    1 + ((value * 2 + index) % 2),
                    1 + ((value + 3) % 4),
                    1 + ((index + value) % 2),
                    1 + ((value * 3 + index) % 3),
                    1 + ((value + index * 2 + 1) % 2)
                ].forEach(append)
                append(1)
            }
        }

        [1, 1, 2, 1, 1, 1, 2].forEach(append)
        return segments
    }

    private var ticketBackground: some View {
        ZStack {
            style.paper

            LinearGradient(
                colors: [
                    Color.white.opacity(0.35),
                    style.accent.opacity(0.045),
                    Color.black.opacity(0.025)
                ],
                startPoint: .topLeading,
                endPoint: .bottomTrailing
            )

            if !model.compactFields.stampAssetName.isEmpty {
                Image(model.compactFields.stampAssetName)
                    .resizable()
                    .scaledToFit()
                    .opacity(model.compactFields.stampOpacity)
                    .rotationEffect(.degrees(model.compactFields.stampRotationDegrees))
                    .frame(width: 86, height: 86)
                    .frame(maxWidth: .infinity, maxHeight: .infinity, alignment: .bottomTrailing)
                    .padding(.trailing, 10)
                    .padding(.bottom, 10)
                    .accessibilityHidden(true)
            }
        }
    }

    private var divider: some View {
        Rectangle()
            .fill(style.separator)
            .frame(height: 1)
    }

    private var languageUsesCJK: Bool {
        switch model.language {
        case .zhHans, .zhHant, .ja, .ko:
            return true
        case .en, .de, .es, .pt, .fr, .ru, .it:
            return false
        }
    }

    @ViewBuilder
    private func optionalText(_ value: String) -> some View {
        if !value.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            Text(value)
        }
    }

    private func moreChipText(count: Int) -> String? {
        let format = model.labels.moreFormat.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !format.isEmpty else { return nil }
        return String(format: format, locale: Locale(identifier: model.language.rawValue), arguments: [count])
    }
}
