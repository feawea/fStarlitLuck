import SwiftUI

struct ZodiacTicketProjectRootView: View {
    var body: some View {
        ZStack {
            Color(red: 0.92, green: 0.93, blue: 0.96)
                .ignoresSafeArea()

            ZodiacTicketView(model: ZodiacTicketProjectSeed.emptyModel)
                .padding(24)
        }
    }
}

enum ZodiacTicketProjectSeed {
    static let emptyLabels = ZodiacTicketLabels(
        title: "",
        serialNumber: "",
        inspiration: "",
        good: "",
        caution: "",
        overview: "",
        triggerFocus: "",
        triggerBadge: "",
        moon: "",
        asc: "",
        moreFormat: ""
    )

    static let emptyFields = ZodiacTicketCompactFields(
        sunMoonAscLine: "",
        inspirationShort: "",
        luckyTitleShort: "",
        luckyOneLineShort: "",
        luckyMoreCount: 0,
        cautionTitleShort: "",
        cautionOneLineShort: "",
        cautionMoreCount: 0,
        focusSummaryShort: "",
        curiosityQuestionShort: "",
        plainReasonShort: "",
        triggerHookLine: "",
        revealLine1: "",
        revealLine2: nil,
        didTearTodayKey: "",
        stampAssetName: "",
        stampOpacity: 0,
        stampRotationDegrees: 0
    )

    static let emptyModel = ZodiacTicketViewModel(
        language: .en,
        displayName: "",
        labels: emptyLabels,
        bigThree: nil,
        compactFields: emptyFields,
        serialNumber: "",
        tearCount: 0
    )
}

