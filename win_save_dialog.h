#pragma once
#include <string>
#include <vector>
#include <windows.h>
#include <commdlg.h>

// --- UTF-8 <-> UTF-16 helpers ---
inline std::wstring utf8_to_wide(const std::string& s) {
    if (s.empty()) return L"";
    int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), nullptr, 0);
    std::wstring w(n, 0);
    MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), &w[0], n);
    return w;
}
inline std::string wide_to_utf8(const std::wstring& w) {
    if (w.empty()) return "";
    int n = WideCharToMultiByte(CP_UTF8, 0, w.c_str(), (int)w.size(), nullptr, 0, nullptr, nullptr);
    std::string s(n, 0);
    WideCharToMultiByte(CP_UTF8, 0, w.c_str(), (int)w.size(), &s[0], n, nullptr, nullptr);
    return s;
}

// فیلترها را به فرم موردنیاز Win32 می‌سازد:  "JSON (*.json)\0*.json\0All (*.*)\0*.*\0\0"
inline std::wstring build_filter(const std::vector<std::pair<std::wstring,std::wstring>>& specs){
    std::wstring f;
    for (auto& sp : specs) {
        f.append(sp.first);  f.push_back(L'\0');
        f.append(sp.second); f.push_back(L'\0');
    }
    f.push_back(L'\0');
    return f;
}

// title: عنوان پنجره (UTF-8)
// suggested_name: نام پیشنهادی فایل (مثلاً "circuit.json")
// default_ext: پسوند پیش‌فرض بدون نقطه (مثلاً "json")
inline std::string show_save_dialog_win(const std::string& title = "Save Circuit",
                                        const std::string& suggested_name = "circuit.json",
                                        const std::string& default_ext = "json")
{
    // بافر مسیر خروجی (UTF-16)
    std::wstring fileBuf(1024, L'\0');
    auto wTitle = utf8_to_wide(title);
    auto wDefExt = utf8_to_wide(default_ext);
    auto wSuggested = utf8_to_wide(suggested_name);

    // اگر نام پیشنهادی داری، داخل بافر بذار
    if (!wSuggested.empty()) {
        size_t n = std::min(wSuggested.size(), fileBuf.size()-1);
        std::copy_n(wSuggested.c_str(), n, &fileBuf[0]);
        fileBuf[n] = L'\0';
    }

    // فیلترها: JSON و All files
    auto filter = build_filter({
                                       {L"JSON Files (*.json)", L"*.json"},
                                       {L"All Files (*.*)",     L"*.*"}
                               });

    OPENFILENAMEW ofn; ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize  = sizeof(ofn);
    ofn.hwndOwner    = nullptr;               // اگر HWND داری می‌تونی اینجا بدهی
    ofn.lpstrFilter  = filter.c_str();
    ofn.lpstrFile    = &fileBuf[0];
    ofn.nMaxFile     = (DWORD)fileBuf.size();
    ofn.lpstrTitle   = wTitle.empty()? nullptr : wTitle.c_str();
    ofn.Flags        = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;
    ofn.lpstrDefExt  = wDefExt.empty()? nullptr : wDefExt.c_str();

    if (GetSaveFileNameW(&ofn)) {
        return wide_to_utf8(std::wstring(ofn.lpstrFile));
    }
    return std::string(); // خالی = کاربر Cancel زده یا خطا
}
