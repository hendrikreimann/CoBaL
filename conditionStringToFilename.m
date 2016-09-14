function clarified_string = conditionStringToFilename(condition_string)
    proxy_string = condition_string;
    proxy_string = strrep(proxy_string, 'LEFT', 'stanceL');
    proxy_string = strrep(proxy_string, 'RIGHT', 'stanceR');
    proxy_string = strrep(proxy_string, 'POSITIVE', 'illuL');
    proxy_string = strrep(proxy_string, 'NEGATIVE', 'illuR');
    proxy_string = strrep(proxy_string, 'ONE', 'step1');
    proxy_string = strrep(proxy_string, 'TWO', 'step2');
    clarified_string = proxy_string;
end