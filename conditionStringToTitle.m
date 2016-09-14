function clarified_string = conditionStringToTitle(condition_string)
    proxy_string = condition_string;
    proxy_string = strrep(proxy_string, 'LEFT', 'stance foot LEFT');
    proxy_string = strrep(proxy_string, 'RIGHT', 'stance foot RIGHT');
    proxy_string = strrep(proxy_string, 'POSITIVE', 'illusion LEFT');
    proxy_string = strrep(proxy_string, 'NEGATIVE', 'illusion RIGHT');
    clarified_string = proxy_string;
end