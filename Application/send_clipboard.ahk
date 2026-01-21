#Requires AutoHotkey v2.0+

; Send content of CLIPBOARD to destination as keystrokes

; CTRL + C to copy text to be sent to destination

; CTRL + SHIFT + V to send contents of clipboard to destination

^+v::{
    Sleep 2000 

    varText := A_Clipboard
    varText := StrReplace(varText, "`r`n", "`n")  ; Normalize line endings

    lines := StrSplit(varText, "`n")
    cleanText := ""
    for line in lines {
        cleanText .= Trim(line) . "`n"
    }

    ; Remove the final newline
    cleanText := SubStr(cleanText, 1, -1)

    ; Send each character one by one
    for char in StrSplit(cleanText) {
        if (char = "`n") {
            SendInput("{Enter}")
        } else {
            SendInput("{Text}" . char)
        }
        Sleep(30)
    }
}