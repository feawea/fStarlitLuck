# ZodiacTicketLocalProject

此目錄是「繁星籤」 的 Xcode 專案備份版，不包含詳情文本、不包含現有票文文本內容。

## 授權與致謝

本備份專案採用 GNU Affero General Public License v3.0 授權，完整英文條款見根目錄 `LICENSE`。

本專案使用 Swiss Ephemeris 的 Swift bridge、C source、headers 與 `.se1` ephemeris data。Swiss Ephemeris 由 Astrodienst AG 發布，作者為 Dieter Koch 與 Alois Treindl；使用、再發布或商業封閉源使用時，請依 Swiss Ephemeris 官方授權條款與 AGPL v3 要求處理。

## 範圍

- `ZodiacTicketLocalProject.xcodeproj`：可單獨打開的 iOS SwiftUI project。
- `ZodiacTicketLocalProject/ZodiacTicketModels.swift`：ticket 所需的占星輸入、Big Three、票面 compact fields 與 view model。
- `ZodiacTicketLocalProject/ZodiacTicketLocalization.swift`：本地 JSON / bundle JSON 載入、模板渲染與 labels 組裝。
- `ZodiacTicketLocalProject/ZodiacTicketComposer.swift`：從輸入資料組出 compact fields / view model。
- `ZodiacTicketLocalProject/ZodiacTicketView.swift`：只呈現 ticket 本體的 SwiftUI view。
- `ZodiacTicketLocalProject/ZodiacTicketLocalProjectApp.swift`：獨立 project app 入口。
- `ZodiacTicketLocalProject/SwissEphemerisBridge.swift`：Swiss Ephemeris Swift service wrapper。
- `ZodiacTicketLocalProject/SwissEphAdapter.swift` 與 `SwissEphAdapterImpl.swift`：星座籤計算層可接入的 adapter。
- `Vendor/SwissEphemeris/src/`：Swiss Ephemeris C source 與 headers。
- `Vendor/SwissEphemeris/ephe/`：Swiss Ephemeris `.se1` data。
- `Resources/zodiac_ticket_text.empty.json`：空白本地化資源模板。
- `Resources/zodiac_ticket_keys.json`：需要由外部資源提供的 key 清單。

## 不包含

- 不包含 `Resources/Zodiac/zodiac_ticket_text.json` 的實際文案。
- 不包含 `Resources/Zodiac/Detail` 與任何詳情頁文案。
- 不包含測試 target。

## 接入方式

打開 `ZodiacTicketLocalProject.xcodeproj` 後可作為獨立本地 project 查看或接入。若要顯示文字，需由外部提供符合 `zodiac_ticket_text.empty.json` 結構的本地化 JSON；若不提供文字資源，composer 與 view 仍可建立，但文字欄位會保持空值。
