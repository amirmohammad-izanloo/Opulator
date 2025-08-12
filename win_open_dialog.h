#pragma once
#include <string>
#include <vector>
#include <windows.h>
#include <commdlg.h>

inline std::wstring build_filter(const std::wstring& title, const std::wstring& pattern){
    // "JSON Files (*.json)\0*.json\0\0"
    std::wstring f = title; f.push_back(L'\0');
    f += pattern; f.push_back(L'\0'); f.push_back(L'\0');
    return f;
}

inline std::string show_open_dialog_win(const std::string& title="Open Circuit",
                                        const std::string& defExt="json")
{
    std::wstring fileBuf(1024, L'\0');
    auto wTitle = utf8_to_wide(title);
    auto filter = build_filter(L"JSON Files (*.json)", L"*.json");

    OPENFILENAMEW ofn; ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize  = sizeof(ofn);
    ofn.hwndOwner    = nullptr;
    ofn.lpstrFilter  = filter.c_str();
    ofn.lpstrFile    = &fileBuf[0];
    ofn.nMaxFile     = (DWORD)fileBuf.size();
    ofn.lpstrTitle   = wTitle.empty()? nullptr : wTitle.c_str();
    ofn.Flags        = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST;
    ofn.lpstrDefExt  = utf8_to_wide(defExt).c_str();

    if (GetOpenFileNameW(&ofn)) {
        return wide_to_utf8(std::wstring(ofn.lpstrFile));
    }
    return {};
}
